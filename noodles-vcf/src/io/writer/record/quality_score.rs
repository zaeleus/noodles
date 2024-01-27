use std::io::{self, Write};

use super::MISSING;

pub(super) fn write_quality_score<W>(writer: &mut W, quality_score: Option<f32>) -> io::Result<()>
where
    W: Write,
{
    if let Some(n) = quality_score {
        write!(writer, "{n}")?;
    } else {
        writer.write_all(MISSING)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_quality_score() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, quality_score: Option<f32>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_quality_score(buf, quality_score)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, b".")?;
        t(&mut buf, Some(0.0), b"0")?;
        t(&mut buf, Some(8.13), b"8.13")?;

        Ok(())
    }
}
