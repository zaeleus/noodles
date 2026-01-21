mod array;
mod character;
mod float;
mod integer;
mod string;

use std::{
    error, fmt,
    io::{self, Write},
};

use self::{
    array::write_array, character::write_character, float::write_float, integer::write_integer,
    string::write_string,
};
use crate::variant::record::info::field::Value;

/// An error returns when an info field value fails to write.
#[derive(Debug)]
pub enum WriteError {
    // I/O error.
    Io(io::Error),
    /// The integer is invalid.
    InvalidInteger(integer::WriteError),
    /// The array is invalid.
    InvalidArray(array::WriteError),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidInteger(e) => Some(e),
            Self::InvalidArray(e) => Some(e),
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidInteger(_) => write!(f, "invalid integer"),
            Self::InvalidArray(_) => write!(f, "invalid array"),
        }
    }
}

pub(super) fn write_value<W>(writer: &mut W, value: &Value) -> Result<(), WriteError>
where
    W: Write,
{
    match value {
        Value::Integer(n) => write_integer(writer, *n).map_err(WriteError::InvalidInteger),
        Value::Float(n) => write_float(writer, *n).map_err(WriteError::Io),
        Value::Flag => Ok(()),
        Value::Character(c) => write_character(writer, *c).map_err(WriteError::Io),
        Value::String(s) => write_string(writer, s).map_err(WriteError::Io),
        Value::Array(array) => write_array(writer, array).map_err(WriteError::InvalidArray),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_value() -> Result<(), WriteError> {
        use crate::variant::record_buf::info::field::Value as ValueBuf;

        fn t(buf: &mut Vec<u8>, value: &ValueBuf, expected: &[u8]) -> Result<(), WriteError> {
            buf.clear();
            write_value(buf, &Value::from(value))?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &ValueBuf::from(0), b"0")?;
        t(&mut buf, &ValueBuf::from(0.0), b"0")?;
        t(&mut buf, &ValueBuf::Flag, b"")?;
        t(&mut buf, &ValueBuf::from('n'), b"n")?;
        t(&mut buf, &ValueBuf::from("noodles"), b"noodles")?;
        t(&mut buf, &ValueBuf::from(vec![Some(8)]), b"8")?;

        Ok(())
    }
}
