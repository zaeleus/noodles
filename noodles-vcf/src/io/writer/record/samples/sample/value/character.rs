use std::io::{self, Write};

pub(super) fn write_character<W>(writer: &mut W, c: char) -> io::Result<()>
where
    W: Write,
{
    write!(writer, "{c}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_character() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_character(&mut buf, 'n')?;
        assert_eq!(buf, b"n");

        buf.clear();
        write_character(&mut buf, ';')?;
        assert_eq!(buf, b";"); // FIXME

        Ok(())
    }
}
