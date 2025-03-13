use std::io::{self, Write};

use super::write_missing;

pub(super) fn write_score<W>(writer: &mut W, score: Option<f32>) -> io::Result<()>
where
    W: Write,
{
    if let Some(n) = score {
        write!(writer, "{n}")
    } else {
        write_missing(writer)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_score() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, score: Option<f32>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_score(buf, score)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, b".")?;
        t(&mut buf, Some(0.0), b"0")?;

        Ok(())
    }
}
