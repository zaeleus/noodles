use std::{io, mem};

use bytes::Buf;
use noodles_sam::{
    self as sam,
    record::cigar::{op::Kind, Op},
};

pub(super) fn get_cigar<B>(
    buf: &mut B,
    cigar: &mut sam::record::Cigar,
    n_cigar_op: usize,
) -> io::Result<()>
where
    B: Buf,
{
    if buf.remaining() < mem::size_of::<u32>() * n_cigar_op {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    cigar.clear();

    for _ in 0..n_cigar_op {
        let op = decode_op(buf.get_u32_le())?;
        cigar.as_mut().push(op);
    }

    Ok(())
}

fn decode_op(n: u32) -> io::Result<Op> {
    let kind = match n & 0x0f {
        0 => Kind::Match,
        1 => Kind::Insertion,
        2 => Kind::Deletion,
        3 => Kind::Skip,
        4 => Kind::SoftClip,
        5 => Kind::HardClip,
        6 => Kind::Pad,
        7 => Kind::SequenceMatch,
        8 => Kind::SequenceMismatch,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid CIGAR op kind",
            ))
        }
    };

    let len = usize::try_from(n >> 4).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(Op::new(kind, len))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_cigar() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            mut src: &[u8],
            actual: &mut sam::record::Cigar,
            n_cigar_op: usize,
            expected: &sam::record::Cigar,
        ) -> io::Result<()> {
            get_cigar(&mut src, actual, n_cigar_op)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        let mut buf = sam::record::Cigar::default();

        t(&[], &mut buf, 0, &sam::record::Cigar::default())?;
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
        assert_eq!(decode_op(0x11)?, Op::new(Kind::Insertion, 1));
        assert_eq!(decode_op(0x12)?, Op::new(Kind::Deletion, 1));
        assert_eq!(decode_op(0x13)?, Op::new(Kind::Skip, 1));
        assert_eq!(decode_op(0x14)?, Op::new(Kind::SoftClip, 1));
        assert_eq!(decode_op(0x15)?, Op::new(Kind::HardClip, 1));
        assert_eq!(decode_op(0x16)?, Op::new(Kind::Pad, 1));
        assert_eq!(decode_op(0x17)?, Op::new(Kind::SequenceMatch, 1));
        assert_eq!(decode_op(0x18)?, Op::new(Kind::SequenceMismatch, 1));

        assert!(matches!(
            decode_op(0x19),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
