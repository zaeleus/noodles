use std::io::{self, Write};

use noodles_core::Position;

use crate::io::writer::num;

pub(super) fn write_position<W>(writer: &mut W, position: Position) -> io::Result<()>
where
    W: Write,
{
    let n = usize::from(position);
    num::write_usize(writer, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_position() -> io::Result<()> {
        let mut buf = Vec::new();
        write_position(&mut buf, Position::MIN)?;
        assert_eq!(buf, b"1");
        Ok(())
    }
}
