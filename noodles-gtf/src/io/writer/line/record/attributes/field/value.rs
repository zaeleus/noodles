use std::io::{self, Write};

pub(super) fn write_value<W>(writer: &mut W, value: &str) -> io::Result<()>
where
    W: Write,
{
    write!(writer, r#""{}""#, value)
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

        Ok(())
    }
}
