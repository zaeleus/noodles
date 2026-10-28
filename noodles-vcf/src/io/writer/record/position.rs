use std::io::{self, Write};

use noodles_core::Position;

pub(super) fn write_position<W>(writer: &mut W, position: Option<Position>) -> io::Result<()>
where
    W: Write,
{
    let n = position.map(usize::from).unwrap_or_default();

    // § 1.6.1.2 "POS" (2026-02-25): "Integer, Required" and § 1.3 "Data types" (2026-02-25):
    // "Integer (32-bit, signed)..."
    if let Ok(m) = i32::try_from(n) {
        write!(writer, "{m}")
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid position: expected <= {}, got {n}", i32::MAX),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_position() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_position(&mut buf, None)?;
        assert_eq!(buf, b"0");

        buf.clear();
        write_position(&mut buf, Some(Position::MIN))?;
        assert_eq!(buf, b"1");

        #[cfg(not(target_pointer_width = "16"))]
        {
            buf.clear();
            write_position(&mut buf, const { Position::new(i32::MAX as usize) })?;
            assert_eq!(buf, b"2147483647");
        }

        #[cfg(not(any(target_pointer_width = "16", target_pointer_width = "32")))]
        {
            buf.clear();
            assert!(matches!(
                write_position(&mut buf, const { Position::new(1 << 31) }),
                Err(e) if e.kind() == io::ErrorKind::InvalidInput
            ));

            buf.clear();
            assert!(matches!(
                write_position(&mut buf, Some(Position::MAX)),
                Err(e) if e.kind() == io::ErrorKind::InvalidInput
            ));
        }

        Ok(())
    }
}
