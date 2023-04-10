mod kind;

use std::io;

use noodles_sam::record::cigar::Op;

use self::kind::decode_kind;

pub(crate) fn decode_op(n: u32) -> io::Result<Op> {
    let kind = decode_kind(n)?;
    let len = usize::try_from(n >> 4).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    Ok(Op::new(kind, len))
}

#[cfg(test)]
mod tests {
    use noodles_sam::record::cigar::op::Kind;

    use super::*;

    #[test]
    fn test_decode_op() -> io::Result<()> {
        assert_eq!(decode_op(0x10)?, Op::new(Kind::Match, 1));

        assert!(matches!(
            decode_op(0x19),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
