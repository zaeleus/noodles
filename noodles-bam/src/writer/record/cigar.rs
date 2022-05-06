use std::io;

use bytes::BufMut;
use noodles_sam::record::{
    cigar::{op::Kind, Op},
    Cigar,
};

pub fn put_cigar<B>(dst: &mut B, cigar: &Cigar) -> io::Result<()>
where
    B: BufMut,
{
    for &op in cigar.as_ref() {
        let n = encode_op(op)?;
        dst.put_u32_le(n);
    }

    Ok(())
}

fn encode_op(op: Op) -> io::Result<u32> {
    const MAX_LENGTH: u32 = (1 << 28) - 1;

    let len =
        u32::try_from(op.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    if len <= MAX_LENGTH {
        let k = encode_kind(op.kind());
        Ok(len << 4 | k)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid CIGAR op length",
        ))
    }
}

fn encode_kind(kind: Kind) -> u32 {
    match kind {
        Kind::Match => 0,
        Kind::Insertion => 1,
        Kind::Deletion => 2,
        Kind::Skip => 3,
        Kind::SoftClip => 4,
        Kind::HardClip => 5,
        Kind::Pad => 6,
        Kind::SequenceMatch => 7,
        Kind::SequenceMismatch => 8,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_cigar() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, cigar: &Cigar, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            put_cigar(buf, cigar)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Cigar::default(), &[])?;
        t(&mut buf, &"4M".parse()?, &[0x40, 0x00, 0x00, 0x00])?;
        t(
            &mut buf,
            &"4M2H".parse()?,
            &[0x40, 0x00, 0x00, 0x00, 0x25, 0x00, 0x00, 0x00],
        )?;

        Ok(())
    }

    #[test]
    fn test_encode_op() -> io::Result<()> {
        let op = Op::new(Kind::Match, 1);
        assert_eq!(encode_op(op)?, 0x10);

        let op = Op::new(Kind::Match, 1 << 28);
        assert!(matches!(
            encode_op(op),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }

    #[test]
    fn test_encode_kind() {
        assert_eq!(encode_kind(Kind::Match), 0);
        assert_eq!(encode_kind(Kind::Insertion), 1);
        assert_eq!(encode_kind(Kind::Deletion), 2);
        assert_eq!(encode_kind(Kind::Skip), 3);
        assert_eq!(encode_kind(Kind::SoftClip), 4);
        assert_eq!(encode_kind(Kind::HardClip), 5);
        assert_eq!(encode_kind(Kind::Pad), 6);
        assert_eq!(encode_kind(Kind::SequenceMatch), 7);
        assert_eq!(encode_kind(Kind::SequenceMismatch), 8);
    }
}
