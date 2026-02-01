mod value;

use std::{
    error, fmt,
    io::{self, Write},
};

use self::value::write_value;
use crate::{Header, io::writer::record::MISSING, variant::record::samples::Sample};

/// An error returns when a record sample fail to write.
#[derive(Debug)]
pub enum WriteError {
    // I/O error.
    Io(io::Error),
    // A value is invalid.
    InvalidValue(value::WriteError),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidValue(e) => Some(e),
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidValue(_) => write!(f, "invalid value"),
        }
    }
}

pub(super) fn write_sample<W, S>(
    writer: &mut W,
    header: &Header,
    sample: S,
) -> Result<(), WriteError>
where
    W: Write,
    S: Sample,
{
    const DELIMITER: &[u8] = b":";

    for (i, result) in sample.iter(header).enumerate() {
        let (_, value) = result.map_err(WriteError::Io)?;

        if i > 0 {
            writer.write_all(DELIMITER).map_err(WriteError::Io)?;
        }

        match value {
            Some(v) => write_value(writer, header, &v).map_err(WriteError::InvalidValue)?,
            None => writer.write_all(MISSING).map_err(WriteError::Io)?,
        }
    }

    Ok(())
}
