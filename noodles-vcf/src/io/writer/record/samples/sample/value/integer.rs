use std::io::{self, Write};

pub(super) fn write_integer<W>(writer: &mut W, n: i32) -> io::Result<()>
where
    W: Write,
{
    write!(writer, "{n}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_integer() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, n: i32, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_integer(buf, n)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, i32::MIN, b"-2147483648")?;
        t(&mut buf, 0, b"0")?;
        t(&mut buf, i32::MAX, b"2147483647")?;

        Ok(())
    }
}
