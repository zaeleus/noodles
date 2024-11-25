use std::io::{self, Write};

use noodles_core::Position;

pub(super) fn write_position<W>(writer: &mut W, position: Position) -> io::Result<()>
where
    W: Write,
{
    write!(writer, "{}", position)
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
