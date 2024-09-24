use std::io::{self, Write};

pub(super) fn write_integer<W>(writer: &mut W, n: i32) -> io::Result<()>
where
    W: Write,
{
    if is_valid(n) {
        write!(writer, "{n}")
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid integer",
        ))
    }
}

// § 1.3 "Data types" (2024-06-28): "For the Integer type, the values from -2^31 to -2^31 + 7
// cannot be stored in the binary version and therefore are disallowed in both VCF and BCF..."
fn is_valid(n: i32) -> bool {
    const MIN: i32 = i32::MIN + 7;
    n > MIN
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

        t(&mut buf, i32::MIN + 8, b"-2147483640")?;
        t(&mut buf, 0, b"0")?;
        t(&mut buf, i32::MAX, b"2147483647")?;

        buf.clear();
        assert!(matches!(
            write_integer(&mut buf, i32::MIN),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
