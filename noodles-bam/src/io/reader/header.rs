//! BAM header reader.

pub(crate) mod magic_number;
mod reference_sequences;
pub mod sam_header;

use std::io::{self, BufRead, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_sam::{self as sam, header::ReferenceSequences};

use self::{magic_number::read_magic_number, reference_sequences::read_reference_sequences};
use crate::MAGIC_NUMBER;

/// A BAM header reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: Read,
{
    pub(super) fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads the magic number.
    pub fn read_magic_number(&mut self) -> io::Result<[u8; MAGIC_NUMBER.len()]> {
        read_magic_number(&mut self.inner)
    }

    /// Returns a SAM header reader.
    ///
    /// The caller is responsible of discarding any extra padding in the header text, e.g., using
    /// [`sam_header::Reader::discard_to_end`].
    pub fn raw_sam_header_reader(&mut self) -> io::Result<sam_header::Reader<&mut R>> {
        let len = self.inner.read_u32::<LittleEndian>().map(u64::from)?;
        Ok(sam_header::Reader::new(&mut self.inner, len))
    }

    /// Reads the reference sequences.
    pub fn read_reference_sequences(&mut self) -> io::Result<ReferenceSequences> {
        read_reference_sequences(&mut self.inner)
    }
}

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<sam::Header>
where
    R: Read,
{
    let mut header_reader = Reader::new(reader);
    read_header_inner(&mut header_reader)
}

fn read_header_inner<R>(reader: &mut Reader<R>) -> io::Result<sam::Header>
where
    R: Read,
{
    reader
        .read_magic_number()
        .and_then(magic_number::validate)?;

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

    fn put_u32_le(dst: &mut Vec<u8>, n: u32) {
        dst.extend(n.to_le_bytes());
    }

    #[test]
    fn test_read_header() -> io::Result<()> {
        let mut src = Vec::new();
        src.extend(MAGIC_NUMBER); // magic
        put_u32_le(&mut src, 27); // l_text
        src.extend(b"@HD\tVN:1.6\n@SQ\tSN:sq0\tLN:8\n"); // text
        put_u32_le(&mut src, 1); // n_ref
        put_u32_le(&mut src, 4); // ref[0].l_name
        src.extend(b"sq0\x00"); // ref[0].name
        put_u32_le(&mut src, 8); // ref[0].l_ref

        let mut reader = &src[..];
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
        let mut src = Vec::new();
        src.extend(MAGIC_NUMBER); // magic
        put_u32_le(&mut src, 11); // l_text
        src.extend(b"@HD\tVN:1.6\n"); // text
        put_u32_le(&mut src, 1); // n_ref
        put_u32_le(&mut src, 4); // ref[0].l_name
        src.extend(b"sq0\x00"); // ref[0].name
        put_u32_le(&mut src, 8); // ref[0].l_ref

        let mut reader = &src[..];
        let actual = read_header(&mut reader)?;

        let expected = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence("sq0", Map::<map::ReferenceSequence>::new(SQ0_LN))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
