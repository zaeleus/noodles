use std::io::{self, Write};

pub(super) fn write_float<W>(writer: &mut W, n: f32) -> io::Result<()>
where
    W: Write,
{
    write!(writer, "{n}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_float() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, n: f32, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_float(buf, n)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 0.0, b"0")?;
        t(&mut buf, 0.333, b"0.333")?;
        t(&mut buf, f32::NAN, b"NaN")?;

        Ok(())
    }
}
