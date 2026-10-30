use std::{
    error, fmt,
    io::{self, Write},
};

use noodles_core::Position;

/// An error returned when a position fails to write.
#[derive(Debug)]
pub enum WriteError {
    /// I/O error.
    Io(io::Error),
    /// The position is invalid.
    InvalidPosition(Option<Position>),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidPosition(_) => None,
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidPosition(position) => {
                let n = position.map(usize::from).unwrap_or_default();
                write!(f, "invalid position: expected <= {}, got {n}", i32::MAX)
            }
        }
    }
}

pub(super) fn write_position<W>(
    writer: &mut W,
    position: Option<Position>,
) -> Result<(), WriteError>
where
    W: Write,
{
    let n = position.map(usize::from).unwrap_or_default();

    // § 1.6.1.2 "POS" (2026-02-25): "Integer, Required" and § 1.3 "Data types" (2026-02-25):
    // "Integer (32-bit, signed)..."
    if let Ok(m) = i32::try_from(n) {
        write!(writer, "{m}").map_err(WriteError::Io)
    } else {
        Err(WriteError::InvalidPosition(position))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_position() -> Result<(), WriteError> {
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
                Err(WriteError::InvalidPosition(_))
            ));

            buf.clear();
            assert!(matches!(
                write_position(&mut buf, Some(Position::MAX)),
                Err(WriteError::InvalidPosition(_))
            ));
        }

        Ok(())
    }
}
