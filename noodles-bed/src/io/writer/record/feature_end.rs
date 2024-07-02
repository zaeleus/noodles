use std::io::{self, Write};

use lexical_core::FormattedSize;
use noodles_core::Position;

pub(super) fn write_feature_end<W>(writer: &mut W, position: Option<Position>) -> io::Result<()>
where
    W: Write,
{
    const MISSING: &[u8] = b"0";

    if let Some(position) = position {
        let n = usize::from(position);
        let mut dst = [0; usize::FORMATTED_SIZE_DECIMAL];
        let buf = lexical_core::write(n, &mut dst);
        writer.write_all(buf)
    } else {
        writer.write_all(MISSING)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_feature_start() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_feature_end(&mut buf, None)?;
        assert_eq!(buf, b"0");

        buf.clear();
        write_feature_end(&mut buf, Some(Position::MIN))?;
        assert_eq!(buf, b"1");

        Ok(())
    }
}
