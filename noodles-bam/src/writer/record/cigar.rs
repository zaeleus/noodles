use std::io;

use bytes::BufMut;
use noodles_sam::{self as sam, record::cigar::Op};

pub fn put_cigar<B>(dst: &mut B, cigar: &sam::record::Cigar) -> io::Result<()>
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
        let k = op.kind() as u32;
        Ok(len << 4 | k)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid CIGAR op length",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_cigar() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, cigar: &sam::record::Cigar, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            put_cigar(buf, cigar)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &sam::record::Cigar::default(), &[])?;
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
        use sam::record::cigar::op::Kind;

        let op = Op::new(Kind::Match, 1);
        assert_eq!(encode_op(op)?, 0x10);

        let op = Op::new(Kind::Match, 1 << 28);
        assert!(matches!(
            encode_op(op),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
