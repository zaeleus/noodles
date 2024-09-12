use std::io::{self, Write};

use percent_encoding::percent_encode_byte;

pub(super) fn write_character<W>(writer: &mut W, c: char) -> io::Result<()>
where
    W: Write,
{
    // § 1.2 "Character encoding, non-printable characters and characters with special meaning" (2024-06-28)
    if c.is_ascii_control() || matches!(c, ':' | ';' | '=' | '%' | ',' | '\r' | '\n' | '\t') {
        // SAFETY: `c` is an ASCII character.
        let b = c as u8;
        let s = percent_encode_byte(b);
        writer.write_all(s.as_bytes())
    } else {
        write!(writer, "{c}")
    }
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
        assert_eq!(buf, b"%3B");

        Ok(())
    }
}
