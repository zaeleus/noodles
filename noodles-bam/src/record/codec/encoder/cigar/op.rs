mod kind;

use std::io;

use noodles_sam::alignment::record::cigar::Op;

use self::kind::encode_kind;

pub(super) fn encode_op(op: Op) -> io::Result<u32> {
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

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record::cigar::op::Kind;

    use super::*;

    #[test]
    fn test_encode_op() -> io::Result<()> {
        assert_eq!(encode_op(Op::new(Kind::Match, 1))?, 0x10);

        assert!(matches!(
            encode_op(Op::new(Kind::Match, 1 << 28)),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }
}
