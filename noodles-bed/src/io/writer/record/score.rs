use std::io::{self, Write};

use lexical_core::FormattedSize;

pub(super) fn write_score<W>(writer: &mut W, score: u16) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; u16::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(score, &mut dst);
    writer.write_all(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_score() -> io::Result<()> {
        let mut buf = Vec::new();
        write_score(&mut buf, 0)?;
        assert_eq!(buf, b"0");
        Ok(())
    }
}
