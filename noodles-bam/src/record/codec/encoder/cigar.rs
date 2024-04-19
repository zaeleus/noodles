mod op;

use std::io;

use bytes::BufMut;
use noodles_sam::alignment::record::Cigar;

use self::op::encode_op;

pub fn put_cigar<B, C>(dst: &mut B, cigar: &C) -> io::Result<()>
where
    B: BufMut,
    C: Cigar,
{
    for result in cigar.iter() {
        let op = result?;
        let n = encode_op(op).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        dst.put_u32_le(n);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::{
        record::cigar::{op::Kind, Op},
        record_buf::Cigar as CigarBuf,
    };

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
}
