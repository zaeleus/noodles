mod op;

use std::{io, mem};

use self::op::parse_op;
use crate::record::Cigar;

pub(crate) fn parse_cigar(mut src: &[u8], cigar: &mut Cigar) -> io::Result<()> {
    let mut ops = Vec::from(mem::take(cigar));

    while !src.is_empty() {
        let op = parse_op(&mut src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        ops.push(op);
    }

    *cigar = Cigar::try_from(ops).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cigar() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::cigar::{op::Kind, Op};

        let src = b"1M13N144S";
        let mut cigar = Cigar::default();
        parse_cigar(src, &mut cigar)?;

        let expected = Cigar::try_from(vec![
            Op::new(Kind::Match, 1),
            Op::new(Kind::Skip, 13),
            Op::new(Kind::SoftClip, 144),
        ])?;

        assert_eq!(cigar, expected);

        Ok(())
    }
}
