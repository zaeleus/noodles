use std::{io, mem};

use bytes::Buf;
use noodles_sam::record::{
    cigar::{op::Kind, Op},
    Cigar,
};

pub fn get_cigar<B>(src: &mut B, cigar: &mut Cigar, n_cigar_op: usize) -> io::Result<()>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u32>() * n_cigar_op {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    cigar.clear();

    for _ in 0..n_cigar_op {
        let op = decode_op(src.get_u32_le())?;
        cigar.as_mut().push(op);
    }

    Ok(())
}

fn decode_op(n: u32) -> io::Result<Op> {
    let kind = decode_kind(n)?;
    let len = usize::try_from(n >> 4).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    Ok(Op::new(kind, len))
}

fn decode_kind(n: u32) -> io::Result<Kind> {
    match n & 0x0f {
        0 => Ok(Kind::Match),
        1 => Ok(Kind::Insertion),
        2 => Ok(Kind::Deletion),
        3 => Ok(Kind::Skip),
        4 => Ok(Kind::SoftClip),
        5 => Ok(Kind::HardClip),
        6 => Ok(Kind::Pad),
        7 => Ok(Kind::SequenceMatch),
        8 => Ok(Kind::SequenceMismatch),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid CIGAR op kind",
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_cigar() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            mut src: &[u8],
            actual: &mut Cigar,
            n_cigar_op: usize,
            expected: &Cigar,
        ) -> io::Result<()> {
            get_cigar(&mut src, actual, n_cigar_op)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        let mut buf = Cigar::default();

        t(&[], &mut buf, 0, &Cigar::default())?;
        t(&[0x40, 0x00, 0x00, 0x00], &mut buf, 1, &"4M".parse()?)?;
        t(
            &[0x40, 0x00, 0x00, 0x00, 0x25, 0x00, 0x00, 0x00],
            &mut buf,
            2,
            &"4M2H".parse()?,
        )?;

        Ok(())
    }

    #[test]
    fn test_decode_op() -> io::Result<()> {
        assert_eq!(decode_op(0x10)?, Op::new(Kind::Match, 1));

        assert!(matches!(
            decode_op(0x19),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_decode_kind() -> io::Result<()> {
        assert_eq!(decode_kind(0x00)?, Kind::Match);
        assert_eq!(decode_kind(0x01)?, Kind::Insertion);
        assert_eq!(decode_kind(0x02)?, Kind::Deletion);
        assert_eq!(decode_kind(0x03)?, Kind::Skip);
        assert_eq!(decode_kind(0x04)?, Kind::SoftClip);
        assert_eq!(decode_kind(0x05)?, Kind::HardClip);
        assert_eq!(decode_kind(0x06)?, Kind::Pad);
        assert_eq!(decode_kind(0x07)?, Kind::SequenceMatch);
        assert_eq!(decode_kind(0x08)?, Kind::SequenceMismatch);

        assert!(matches!(
            decode_op(0x09),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
