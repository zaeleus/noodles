mod op;

use std::{io, mem};

use bytes::Buf;
use noodles_sam::record::Cigar;

pub(super) use self::op::decode_op;

pub fn get_cigar<B>(src: &mut B, cigar: &mut Cigar, n_cigar_op: usize) -> io::Result<()>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u32>() * n_cigar_op {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    cigar.clear();

    for _ in 0..n_cigar_op {
        let op = decode_op(src.get_u32_le())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        cigar.as_mut().push(op);
    }

    Ok(())
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
}
