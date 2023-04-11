mod kind;

use std::io;

use self::kind::parse_kind;
use crate::record::cigar::Op;

pub(super) fn parse_op(src: &mut &[u8]) -> io::Result<Op> {
    let len = parse_len(src)?;
    let kind = parse_kind(src)?;
    Ok(Op::new(kind, len))
}

fn parse_len(src: &mut &[u8]) -> io::Result<usize> {
    let (len, i) = lexical_core::parse_partial(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *src = &src[i..];

    Ok(len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_len() {
        assert!(matches!(
            parse_len(&mut &[][..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));
    }
}
