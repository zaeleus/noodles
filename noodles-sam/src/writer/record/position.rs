use std::io::{self, Write};

use noodles_core::Position;

use crate::writer::num;

// § 1.4 "The alignment section: mandatory fields" (2021-06-03): `[0, 2^31-1]`.
const MAX_POSITION_VALUE: usize = (1 << 31) - 1;

pub fn write_position<W>(writer: &mut W, position: Option<Position>) -> io::Result<()>
where
    W: Write,
{
    let pos = position.map(usize::from).unwrap_or_default();

    if pos <= MAX_POSITION_VALUE {
        num::write_usize(writer, pos)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "invalid position: expected position to be <= {}, got {}",
                MAX_POSITION_VALUE, pos
            ),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_position() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, position: Option<Position>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_position(buf, position)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, b"0")?;
        t(&mut buf, Position::new(13), b"13")?;

        buf.clear();
        assert!(matches!(
            write_position(&mut buf, Position::new(MAX_POSITION_VALUE + 1)),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }
}
