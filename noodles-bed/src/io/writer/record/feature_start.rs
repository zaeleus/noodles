use std::io::{self, Write};

use lexical_core::FormattedSize;
use noodles_core::Position;

pub(super) fn write_feature_start<W>(writer: &mut W, position: Position) -> io::Result<()>
where
    W: Write,
{
    // § 1.6.2 "Coordinates: `chromStart`" (2022-01-05): "...`chromStart` must be less than or
    // equal to 2^64 - 1..."
    let n = u64::try_from(usize::from(position))
        .map(|n| n - 1)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let mut dst = [0; u64::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_feature_start() -> io::Result<()> {
        let mut buf = Vec::new();
        write_feature_start(&mut buf, Position::MIN)?;
        assert_eq!(buf, b"0");
        Ok(())
    }
}
