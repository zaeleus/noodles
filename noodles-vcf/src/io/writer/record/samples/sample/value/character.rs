use std::io::{self, Write};

use percent_encoding::percent_encode_byte;

pub(super) fn write_character<W>(writer: &mut W, c: char) -> io::Result<()>
where
    W: Write,
{
    // § 1.2 "Character encoding, non-printable characters and characters with special meaning"
    // (2024-10-09).
    if c.is_ascii_control() || matches!(c, ':' | '%' | ',' | '\r' | '\n' | '\t') {
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
        fn t(buf: &mut Vec<u8>, c: char, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_character(buf, c)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 'n', b"n")?;
        t(&mut buf, ';', b";")?;
        t(&mut buf, '=', b"=")?;
        t(&mut buf, ':', b"%3A")?;

        Ok(())
    }
}
