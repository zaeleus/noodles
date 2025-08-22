use std::io::{self, Write};

use crate::alignment::record::cigar::op::Kind;

pub(super) fn write_kind<W>(writer: &mut W, kind: Kind) -> io::Result<()>
where
    W: Write,
{
    let c = match kind {
        Kind::Match => b'M',
        Kind::Insertion => b'I',
        Kind::Deletion => b'D',
        Kind::Skip => b'N',
        Kind::SoftClip => b'S',
        Kind::HardClip => b'H',
        Kind::Pad => b'P',
        Kind::SequenceMatch => b'=',
        Kind::SequenceMismatch => b'X',
    };

    writer.write_all(&[c])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_kind() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, kind: Kind, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_kind(buf, kind)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, Kind::Match, b"M")?;
        t(&mut buf, Kind::Insertion, b"I")?;
        t(&mut buf, Kind::Deletion, b"D")?;
        t(&mut buf, Kind::Skip, b"N")?;
        t(&mut buf, Kind::SoftClip, b"S")?;
        t(&mut buf, Kind::HardClip, b"H")?;
        t(&mut buf, Kind::Pad, b"P")?;
        t(&mut buf, Kind::SequenceMatch, b"=")?;
        t(&mut buf, Kind::SequenceMismatch, b"X")?;

        Ok(())
    }
}
