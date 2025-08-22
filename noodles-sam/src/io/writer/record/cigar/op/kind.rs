use std::io::{self, Write};

use crate::alignment::record::cigar::op::Kind;

pub(super) fn write_kind<W>(writer: &mut W, kind: Kind) -> io::Result<()>
where
    W: Write,
{
    let c = encode(kind);
    writer.write_all(&[c])
}

fn encode(kind: Kind) -> u8 {
    match kind {
        Kind::Match => b'M',
        Kind::Insertion => b'I',
        Kind::Deletion => b'D',
        Kind::Skip => b'N',
        Kind::SoftClip => b'S',
        Kind::HardClip => b'H',
        Kind::Pad => b'P',
        Kind::SequenceMatch => b'=',
        Kind::SequenceMismatch => b'X',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_kind() -> io::Result<()> {
        let mut buf = Vec::new();
        write_kind(&mut buf, Kind::Match)?;
        assert_eq!(buf, b"M");
        Ok(())
    }

    #[test]
    fn test_encode() {
        assert_eq!(encode(Kind::Match), b'M');
        assert_eq!(encode(Kind::Insertion), b'I');
        assert_eq!(encode(Kind::Deletion), b'D');
        assert_eq!(encode(Kind::Skip), b'N');
        assert_eq!(encode(Kind::SoftClip), b'S');
        assert_eq!(encode(Kind::HardClip), b'H');
        assert_eq!(encode(Kind::Pad), b'P');
        assert_eq!(encode(Kind::SequenceMatch), b'=');
        assert_eq!(encode(Kind::SequenceMismatch), b'X');
    }
}
