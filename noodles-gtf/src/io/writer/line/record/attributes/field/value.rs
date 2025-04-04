use std::io::{self, Write};

use bstr::BStr;

const BACKSLASH: u8 = b'\\';
const QUOTATION_MARK: u8 = b'"';

pub(super) fn write_value<W>(writer: &mut W, value: &BStr) -> io::Result<()>
where
    W: Write,
{
    if requires_escapes(value) {
        write_escaped_string(writer, value)?;
    } else {
        writer.write_all(&[QUOTATION_MARK])?;
        writer.write_all(value)?;
        writer.write_all(&[QUOTATION_MARK])?;
    }

    Ok(())
}

fn requires_escapes(s: &BStr) -> bool {
    s.iter().any(|c| matches!(*c, BACKSLASH | QUOTATION_MARK))
}

fn write_escaped_string<W>(writer: &mut W, s: &BStr) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[QUOTATION_MARK])?;

    for c in s.iter().copied() {
        if matches!(c, BACKSLASH | QUOTATION_MARK) {
            writer.write_all(&[BACKSLASH])?;
        }

        writer.write_all(&[c])?;
    }

    writer.write_all(&[QUOTATION_MARK])?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_value() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, value: &BStr, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, BStr::new(b"ndls"), br#""ndls""#)?;
        t(&mut buf, BStr::new(br#"nd\ls"#), br#""nd\\ls""#)?;
        t(&mut buf, BStr::new(br#"nd"ls""#), br#""nd\"ls\"""#)?;

        Ok(())
    }
}
