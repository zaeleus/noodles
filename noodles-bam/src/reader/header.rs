use std::{
    io::{self, BufRead, BufReader, Read},
    num::NonZeroUsize,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_sam::{
    self as sam,
    header::{
        record::value::{
            map::{self, ReferenceSequence},
            Map,
        },
        ReferenceSequences,
    },
};

use super::bytes_with_nul_to_string;
use crate::MAGIC_NUMBER;

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<sam::Header>
where
    R: Read,
{
    read_magic(reader)?;

    let mut header = read_header_inner(reader)?;
    let reference_sequences = read_reference_sequences(reader)?;

    if header.reference_sequences().is_empty() {
        *header.reference_sequences_mut() = reference_sequences;
    } else if !reference_sequences_eq(header.reference_sequences(), &reference_sequences) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "SAM header and binary reference sequence dictionaries mismatch",
        ));
    }

    Ok(header)
}

fn read_magic<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    let mut magic = [0; 4];
    reader.read_exact(&mut magic)?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BAM header",
        ))
    }
}

fn read_header_inner<R>(reader: &mut R) -> io::Result<sam::Header>
where
    R: Read,
{
    const NUL: u8 = 0x00;

    let l_text = reader.read_u32::<LittleEndian>().map(u64::from)?;

    let mut parser = sam::header::Parser::default();

    let mut header_reader = BufReader::new(reader.take(l_text));
    let mut buf = Vec::new();

    while read_line(&mut header_reader, &mut buf)? != 0 {
        // § 4.2 "The BAM format" (2023-05-24): "Plain header text in SAM; not necessarily
        // `NUL`-terminated".
        if buf == [NUL] {
            break;
        }

        parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    Ok(parser.finish())
}

fn read_line<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    dst.clear();

    match reader.read_until(LINE_FEED, dst)? {
        0 => Ok(0),
        n => {
            if dst.ends_with(&[LINE_FEED]) {
                dst.pop();

                if dst.ends_with(&[CARRIAGE_RETURN]) {
                    dst.pop();
                }
            }

            Ok(n)
        }
    }
}

fn read_reference_sequences<R>(reader: &mut R) -> io::Result<ReferenceSequences>
where
    R: Read,
{
    let n_ref = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut reference_sequences = ReferenceSequences::with_capacity(n_ref);

    for _ in 0..n_ref {
        let (name, reference_sequence) = read_reference_sequence(reader)?;
        reference_sequences.insert(name, reference_sequence);
    }

    Ok(reference_sequences)
}

fn read_reference_sequence<R>(
    reader: &mut R,
) -> io::Result<(map::reference_sequence::Name, Map<ReferenceSequence>)>
where
    R: Read,
{
    let l_name = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut c_name = vec![0; l_name];
    reader.read_exact(&mut c_name)?;

    let name = bytes_with_nul_to_string(&c_name).and_then(|name| {
        name.parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let l_ref = reader.read_u32::<LittleEndian>().and_then(|len| {
        usize::try_from(len)
            .and_then(NonZeroUsize::try_from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let reference_sequence = Map::<ReferenceSequence>::new(l_ref);

    Ok((name, reference_sequence))
}

fn reference_sequences_eq(
    header_reference_sequences: &ReferenceSequences,
    binary_reference_sequences: &ReferenceSequences,
) -> bool {
    header_reference_sequences.len() == binary_reference_sequences.len()
        && header_reference_sequences
            .iter()
            .zip(binary_reference_sequences)
            .all(|((h_name, h_map), (b_name, b_map))| {
                h_name == b_name && h_map.length() == b_map.length()
            })
}

#[cfg(test)]
mod tests {
    use bytes::BufMut;
    use noodles_sam::{
        self as sam,
        header::record::value::{
            map::{self, header::Version},
            Map,
        },
    };

    use super::*;

    #[test]
    fn test_read_magic() -> io::Result<()> {
        let data = b"BAM\x01";
        let mut reader = &data[..];
        assert!(read_magic(&mut reader).is_ok());

        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"MThd";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_read_header() -> Result<(), Box<dyn std::error::Error>> {
        let mut data = Vec::new();
        data.put_slice(MAGIC_NUMBER); // magic
        data.put_u32_le(27); // l_text
        data.put_slice(b"@HD\tVN:1.6\n@SQ\tSN:sq0\tLN:8\n"); // text
        data.put_u32_le(1); // n_ref
        data.put_u32_le(4); // ref[0].l_name
        data.put_slice(b"sq0\x00"); // ref[0].name
        data.put_u32_le(8); // ref[0].l_ref

        let mut reader = &data[..];
        let actual = read_header(&mut reader)?;

        let expected = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence(
                "sq0".parse()?,
                Map::<map::ReferenceSequence>::new(NonZeroUsize::try_from(8)?),
            )
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_with_trailing_nul_in_text() -> Result<(), Box<dyn std::error::Error>> {
        let mut data = Vec::new();
        data.put_slice(MAGIC_NUMBER); // magic
        data.put_u32_le(12); // l_text
        data.put_slice(b"@HD\tVN:1.6\n\x00"); // text
        data.put_u32_le(0); // n_ref

        let mut reader = &data[..];
        let actual = read_header(&mut reader)?;

        let expected = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_with_missing_sam_header_reference_sequence_dictionary(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut data = Vec::new();
        data.put_slice(MAGIC_NUMBER); // magic
        data.put_u32_le(11); // l_text
        data.put_slice(b"@HD\tVN:1.6\n"); // text
        data.put_u32_le(1); // n_ref
        data.put_u32_le(4); // ref[0].l_name
        data.put_slice(b"sq0\x00"); // ref[0].name
        data.put_u32_le(8); // ref[0].l_ref

        let mut reader = &data[..];
        let actual = read_header(&mut reader)?;

        let expected = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence(
                "sq0".parse()?,
                Map::<map::ReferenceSequence>::new(NonZeroUsize::try_from(8)?),
            )
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_reference_sequences() -> Result<(), Box<dyn std::error::Error>> {
        let data = [
            0x01, 0x00, 0x00, 0x00, // n_ref = 1
            0x04, 0x00, 0x00, 0x00, // ref[0].l_name = 4
            0x73, 0x71, 0x30, 0x00, // ref[0].name = "sq0\x00"
            0x08, 0x00, 0x00, 0x00, // ref[0].l_ref = 8
        ];

        let mut reader = &data[..];
        let actual = read_reference_sequences(&mut reader)?;

        let expected: ReferenceSequences = [(
            "sq0".parse()?,
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(8)?),
        )]
        .into_iter()
        .collect();

        assert_eq!(actual, expected);

        Ok(())
    }
}
