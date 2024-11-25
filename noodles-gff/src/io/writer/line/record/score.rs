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
        let mut buf = Vec::new();

        buf.clear();
        write_score(&mut buf, None)?;
        assert_eq!(buf, b".");

        buf.clear();
        write_score(&mut buf, Some(0.0))?;
        assert_eq!(buf, b"0");

        Ok(())
    }
}
