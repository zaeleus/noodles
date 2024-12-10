mod reference_sequences;
mod sam_header;

use std::io::{self, BufRead, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_sam::{self as sam, header::ReferenceSequences};

use self::reference_sequences::read_reference_sequences;
use crate::MAGIC_NUMBER;

struct Reader<'r, R> {
    inner: &'r mut R,
}

impl<'r, R> Reader<'r, R>
where
    R: Read,
{
    fn new(inner: &'r mut R) -> Self {
        Self { inner }
    }

    fn read_magic_number(&mut self) -> io::Result<[u8; MAGIC_NUMBER.len()]> {
        let mut buf = [0; MAGIC_NUMBER.len()];
        self.inner.read_exact(&mut buf)?;
        Ok(buf)
    }

    fn raw_sam_header_reader(&mut self) -> io::Result<sam_header::Reader<R>> {
        let len = self.inner.read_u32::<LittleEndian>().map(u64::from)?;
        Ok(sam_header::Reader::new(self.inner, len))
    }

    fn read_reference_sequences(&mut self) -> io::Result<ReferenceSequences> {
        read_reference_sequences(self.inner)
    }
}

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<sam::Header>
where
    R: Read,
{
    let mut header_reader = Reader::new(reader);
    read_header_inner(&mut header_reader)
}

fn read_magic_number<R>(reader: &mut Reader<R>) -> io::Result<()>
where
    R: Read,
{
    let magic_number = reader.read_magic_number()?;

    if magic_number == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BAM header",
        ))
    }
}

fn read_header_inner<R>(reader: &mut Reader<R>) -> io::Result<sam::Header>
where
    R: Read,
{
    read_magic_number(reader)?;

    let mut raw_sam_header_reader = reader.raw_sam_header_reader()?;
    let mut header = read_sam_header(&mut raw_sam_header_reader)?;

    let reference_sequences = reader.read_reference_sequences()?;

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

fn read_sam_header<R>(reader: &mut sam_header::Reader<R>) -> io::Result<sam::Header>
where
    R: Read,
{
    let mut parser = sam::header::Parser::default();

    let mut buf = Vec::new();

    while read_line(reader, &mut buf)? != 0 {
        parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    reader.discard_to_end()?;

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

pub(crate) fn reference_sequences_eq(
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
    use std::num::NonZeroUsize;

    use bytes::BufMut;
    use noodles_sam::{
        self as sam,
        header::record::value::{
            map::{self, header::Version},
            Map,
        },
    };

    use super::*;

    const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
        Some(length) => length,
        None => unreachable!(),
    };

    #[test]
    fn test_read_magic_number() -> io::Result<()> {
        let mut src = &b"BAM\x01"[..];
        let mut reader = Reader::new(&mut src);
        assert!(read_magic_number(&mut reader).is_ok());

        let mut src = &[][..];
        let mut reader = Reader::new(&mut src);
        assert!(matches!(
            read_magic_number(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let mut src = &b"MThd"[..];
        let mut reader = Reader::new(&mut src);
        assert!(matches!(
            read_magic_number(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_read_header() -> io::Result<()> {
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
            .add_reference_sequence("sq0", Map::<map::ReferenceSequence>::new(SQ0_LN))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_with_missing_sam_header_reference_sequence_dictionary() -> io::Result<()> {
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
            .add_reference_sequence("sq0", Map::<map::ReferenceSequence>::new(SQ0_LN))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
