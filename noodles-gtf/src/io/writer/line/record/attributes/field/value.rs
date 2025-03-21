use std::io::{self, Write};

const BACKSLASH: char = '\\';
const QUOTATION_MARK: char = '"';

pub(super) fn write_value<W>(writer: &mut W, value: &str) -> io::Result<()>
where
    W: Write,
{
    if requires_escapes(value) {
        write_escaped_string(writer, value)
    } else {
        write!(writer, "{}{}{}", QUOTATION_MARK, value, QUOTATION_MARK)
    }
}

fn requires_escapes(s: &str) -> bool {
    s.contains(|c| matches!(c, BACKSLASH | QUOTATION_MARK))
}

fn write_escaped_string<W>(writer: &mut W, s: &str) -> io::Result<()>
where
    W: Write,
{
    write!(writer, "{QUOTATION_MARK}")?;

    for c in s.chars() {
        if matches!(c, BACKSLASH | QUOTATION_MARK) {
            write!(writer, "{BACKSLASH}")?;
        }

        write!(writer, "{c}")?;
    }

    write!(writer, "{QUOTATION_MARK}")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_value() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, value: &str, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, "ndls", br#""ndls""#)?;
        t(&mut buf, r#"nd\ls"#, br#""nd\\ls""#)?;
        t(&mut buf, r#"nd"ls""#, br#""nd\"ls\"""#)?;

        Ok(())
    }
}
