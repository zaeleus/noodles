use std::{error, fmt, num};

use bytes::BufMut;
use noodles_core::Position;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EncodeError {
    Invalid(num::TryFromIntError),
}

impl error::Error for EncodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Invalid(e) => Some(e),
        }
    }
}

impl fmt::Display for EncodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(_) => write!(f, "invalid input"),
        }
    }
}

pub(super) fn put_position<B>(dst: &mut B, position: Option<Position>) -> Result<(), EncodeError>
where
    B: BufMut,
{
    const MISSING: i32 = -1;

    let pos = if let Some(position) = position {
        let n = usize::from(position) - 1;
        i32::try_from(n).map_err(EncodeError::Invalid)?
    } else {
        MISSING
    };

    dst.put_i32_le(pos);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_position() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            buf: &mut Vec<u8>,
            position: Option<Position>,
            expected: &[u8],
        ) -> Result<(), EncodeError> {
            buf.clear();
            put_position(buf, position)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, &[0xff, 0xff, 0xff, 0xff])?;
        t(&mut buf, Some(Position::MIN), &[0x00, 0x00, 0x00, 0x00])?;
        t(
            &mut buf,
            Position::try_from(8).map(Some)?,
            &[0x07, 0x00, 0x00, 0x00],
        )?;

        Ok(())
    }
}
