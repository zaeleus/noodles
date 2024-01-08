use std::io;

use bytes::BufMut;
use noodles_sam::alignment::record::{cigar::op::Kind, Cigar};

pub fn put_cigar<B, C>(dst: &mut B, cigar: &C) -> io::Result<()>
where
    B: BufMut,
    C: Cigar,
{
    for result in cigar.iter() {
        let (kind, len) = result?;
        let n = encode_op(kind, len)?;
        dst.put_u32_le(n);
    }

    Ok(())
}

fn encode_op(kind: Kind, len: usize) -> io::Result<u32> {
    const MAX_LENGTH: u32 = (1 << 28) - 1;

    let len = u32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    if len <= MAX_LENGTH {
        let k = encode_kind(kind);
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
    use noodles_sam::{alignment::record::cigar::Op, record::Cigar as CigarBuf};

    use super::*;

    #[test]
    fn test_put_cigar() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, cigar: &CigarBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            put_cigar(buf, cigar)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &CigarBuf::default(), &[])?;
        t(
            &mut buf,
            &[Op::new(Kind::Match, 4)].into_iter().collect(),
            &[0x40, 0x00, 0x00, 0x00],
        )?;
        t(
            &mut buf,
            &[Op::new(Kind::Match, 4), Op::new(Kind::HardClip, 2)]
                .into_iter()
                .collect(),
            &[0x40, 0x00, 0x00, 0x00, 0x25, 0x00, 0x00, 0x00],
        )?;

        Ok(())
    }

    #[test]
    fn test_encode_op() -> io::Result<()> {
        assert_eq!(encode_op(Kind::Match, 1)?, 0x10);

        assert!(matches!(
            encode_op(Kind::Match, 1 << 28),
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
