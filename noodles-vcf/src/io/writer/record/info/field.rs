mod key;
mod value;

use std::{
    error, fmt,
    io::{self, Write},
};

use self::{key::write_key, value::write_value};
use crate::{io::writer::record::MISSING, variant::record::info::field::Value};

/// An error returns when an info field fails to write.
#[derive(Debug)]
pub enum WriteError {
    // I/O error.
    Io(io::Error),
    /// The key is invalid.
    InvalidKey(key::WriteError),
    /// The value is invalid.
    InvalidValue(value::WriteError),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidKey(e) => Some(e),
            Self::InvalidValue(e) => Some(e),
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidKey(_) => write!(f, "invalid key"),
            Self::InvalidValue(_) => write!(f, "invalid value"),
        }
    }
}

pub(super) fn write_field<W>(
    writer: &mut W,
    key: &str,
    value: Option<&Value>,
) -> Result<(), WriteError>
where
    W: Write,
{
    const SEPARATOR: &[u8] = b"=";

    write_key(writer, key).map_err(WriteError::InvalidKey)?;

    match value {
        Some(Value::Flag) => {}
        Some(v) => {
            writer.write_all(SEPARATOR).map_err(WriteError::Io)?;
            write_value(writer, v).map_err(WriteError::InvalidValue)?;
        }
        None => {
            writer.write_all(SEPARATOR).map_err(WriteError::Io)?;
            writer.write_all(MISSING).map_err(WriteError::Io)?;
        }
    }

    Ok(())
}
