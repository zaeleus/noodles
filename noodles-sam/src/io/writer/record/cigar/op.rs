mod kind;

use std::io::{self, Write};

use self::kind::write_kind;
use crate::{alignment::record::cigar::Op, io::writer::num};

pub(super) fn write_op<W>(writer: &mut W, op: Op) -> io::Result<()>
where
    W: Write,
{
    num::write_usize(writer, op.len())?;
    write_kind(writer, op.kind())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::record::cigar::op::Kind;

    #[test]
    fn test_write_op() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_op(&mut buf, Op::new(Kind::Match, 1))?;
        assert_eq!(buf, b"1M");

        buf.clear();
        write_op(&mut buf, Op::new(Kind::Match, 1 << 28))?;
        assert_eq!(buf, b"268435456M");

        Ok(())
    }
}
