mod op;

use std::io;

use noodles_sam::alignment::record::Cigar;

use self::op::encode_op;
use super::num::write_u32_le;

pub fn write_cigar<C>(dst: &mut Vec<u8>, cigar: &C) -> io::Result<()>
where
    C: Cigar,
{
    for result in cigar.iter() {
        let op = result?;
        let n = encode_op(op).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_u32_le(dst, n);
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
            write_cigar(buf, cigar)?;
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
