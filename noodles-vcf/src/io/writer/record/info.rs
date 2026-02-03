mod field;

use std::{
    error, fmt,
    io::{self, Write},
};

use self::field::write_field;
use super::MISSING;
use crate::{Header, variant::record::Info};

/// An error returns when info fields fail to write.
#[derive(Debug)]
pub enum WriteError {
    /// I/O error.
    Io(io::Error),
    /// A field is invalid.
    InvalidField(field::WriteError),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidField(e) => Some(e),
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidField(_) => write!(f, "invalid field"),
        }
    }
}

pub(super) fn write_info<W, I>(writer: &mut W, header: &Header, info: I) -> Result<(), WriteError>
where
    W: Write,
    I: Info,
{
    const DELIMITER: &[u8] = b";";

    if info.is_empty() {
        writer.write_all(MISSING).map_err(WriteError::Io)?;
    } else {
        for (i, result) in info.iter(header).enumerate() {
            let (key, value) = result.map_err(WriteError::Io)?;

            if i > 0 {
                writer.write_all(DELIMITER).map_err(WriteError::Io)?;
            }

            write_field(writer, key, value.as_ref()).map_err(WriteError::InvalidField)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_info() -> Result<(), WriteError> {
        use crate::variant::{
            record::info::field::key,
            record_buf::{Info as InfoBuf, info::field::Value as ValueBuf},
        };

        fn t(
            buf: &mut Vec<u8>,
            header: &Header,
            info: &InfoBuf,
            expected: &[u8],
        ) -> Result<(), WriteError> {
            buf.clear();
            write_info(buf, header, info)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let header = Header::default();
        let mut buf = Vec::new();

        let info = InfoBuf::default();
        t(&mut buf, &header, &info, b".")?;

        let info = [(
            String::from(key::SAMPLES_WITH_DATA_COUNT),
            Some(ValueBuf::from(2)),
        )]
        .into_iter()
        .collect();
        t(&mut buf, &header, &info, b"NS=2")?;

        let info = [
            (
                String::from(key::SAMPLES_WITH_DATA_COUNT),
                Some(ValueBuf::from(2)),
            ),
            (String::from(key::IS_IN_DB_SNP), Some(ValueBuf::Flag)),
        ]
        .into_iter()
        .collect();

        t(&mut buf, &header, &info, b"NS=2;DB")?;

        Ok(())
    }
}
