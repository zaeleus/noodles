use std::io::{self, Write};

use noodles_core::Position;

pub(super) fn write_position<W>(writer: &mut W, position: Option<Position>) -> io::Result<()>
where
    W: Write,
{
    let n = position.map(usize::from).unwrap_or_default();
    write!(writer, "{}", n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_position() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_position(&mut buf, None)?;
        assert_eq!(buf, b"0");

        buf.clear();
        write_position(&mut buf, Some(Position::MIN))?;
        assert_eq!(buf, b"1");

        Ok(())
    }
}
