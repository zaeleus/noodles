use std::{
    io::{self, Read},
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

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<String>
where
    R: Read,
{
    read_magic(reader)?;
    read_raw_header(reader)
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

fn read_raw_header<R>(reader: &mut R) -> io::Result<String>
where
    R: Read,
{
    let l_text = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut text = vec![0; l_text];
    reader.read_exact(&mut text)?;

    // § 4.2 The BAM format (2021-06-03): "Plain header text in SAM; not necessarily
    // NUL-terminated".
    bytes_with_nul_to_string(&text).or_else(|_| {
        String::from_utf8(text).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

pub(super) fn read_reference_sequences<R>(reader: &mut R) -> io::Result<ReferenceSequences>
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

pub(super) fn read_alignment_header<R>(reader: &mut R) -> io::Result<sam::Header>
where
    R: Read,
{
    let header = read_header(reader).and_then(|s| {
        s.parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    read_reference_sequences(reader)?;

    Ok(header)
}

#[cfg(test)]
mod tests {
    use noodles_sam as sam;

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
    fn test_read_raw_header() -> io::Result<()> {
        let expected = "@HD\tVN:1.6\n";

        let data_len = expected.len() as u32;
        let mut data = data_len.to_le_bytes().to_vec();
        data.extend(expected.as_bytes());

        let mut reader = &data[..];
        let actual = read_raw_header(&mut reader)?;

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

    #[test]
    fn test_read_alignment_header() -> Result<(), Box<dyn std::error::Error>> {
        use bytes::BufMut;
        use sam::header::record::value::{
            map::{self, header::Version},
            Map,
        };

        let mut data = Vec::new();
        data.put_slice(MAGIC_NUMBER); // magic
        data.put_u32_le(11); // l_text
        data.put_slice(b"@HD\tVN:1.6\n"); // text
        data.put_u32_le(1); // n_ref
        data.put_u32_le(4); // ref[0].l_name
        data.put_slice(b"sq0\x00"); // ref[0].name
        data.put_u32_le(8); // ref[0].l_ref

        let mut reader = &data[..];
        let actual = read_alignment_header(&mut reader)?;

        let expected = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
