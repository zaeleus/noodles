use std::{
    error, fmt,
    io::{self, Write},
};

use super::{write_character, write_float, write_integer, write_string};
use crate::{io::writer::record::MISSING, variant::record::samples::series::value::Array};

/// An error returns when a sample array value fails to write.
#[derive(Debug)]
pub enum WriteError {
    /// I/O error.
    Io(io::Error),
    /// An integer is invalid.
    InvalidInteger(super::integer::WriteError),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidInteger(e) => Some(e),
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidInteger(_) => write!(f, "invalid integer"),
        }
    }
}

pub(super) fn write_array<W>(writer: &mut W, array: &Array) -> Result<(), WriteError>
where
    W: Write,
{
    match array {
        Array::Integer(values) => {
            for (i, result) in values.iter().enumerate() {
                if i > 0 {
                    write_separator(writer).map_err(WriteError::Io)?;
                }

                if let Some(n) = result.map_err(WriteError::Io)? {
                    write_integer(writer, n).map_err(WriteError::InvalidInteger)?;
                } else {
                    writer.write_all(MISSING).map_err(WriteError::Io)?;
                }
            }
        }
        Array::Float(values) => {
            for (i, result) in values.iter().enumerate() {
                if i > 0 {
                    write_separator(writer).map_err(WriteError::Io)?;
                }

                if let Some(n) = result.map_err(WriteError::Io)? {
                    write_float(writer, n).map_err(WriteError::Io)?;
                } else {
                    writer.write_all(MISSING).map_err(WriteError::Io)?;
                }
            }
        }
        Array::Character(values) => {
            for (i, result) in values.iter().enumerate() {
                if i > 0 {
                    write_separator(writer).map_err(WriteError::Io)?;
                }

                if let Some(c) = result.map_err(WriteError::Io)? {
                    write_character(writer, c).map_err(WriteError::Io)?;
                } else {
                    writer.write_all(MISSING).map_err(WriteError::Io)?;
                }
            }
        }
        Array::String(values) => {
            for (i, result) in values.iter().enumerate() {
                if i > 0 {
                    write_separator(writer).map_err(WriteError::Io)?;
                }

                if let Some(s) = result.map_err(WriteError::Io)? {
                    write_string(writer, &s).map_err(WriteError::Io)?;
                } else {
                    writer.write_all(MISSING).map_err(WriteError::Io)?;
                }
            }
        }
    }

    Ok(())
}

fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: &[u8] = b",";
    writer.write_all(SEPARATOR)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record_buf::samples::sample::value::Array as ArrayBuf;

    #[test]
    fn test_write_array() -> Result<(), WriteError> {
        fn t(buf: &mut Vec<u8>, array: &ArrayBuf, expected: &[u8]) -> Result<(), WriteError> {
            buf.clear();
            write_array(buf, &array.into())?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let array = ArrayBuf::Integer(vec![Some(8)]);
        t(&mut buf, &array, b"8")?;

        let array = ArrayBuf::Integer(vec![Some(8), Some(13), None]);
        t(&mut buf, &array, b"8,13,.")?;

        let array = ArrayBuf::Float(vec![Some(0.333)]);
        t(&mut buf, &array, b"0.333")?;

        let array = ArrayBuf::Float(vec![Some(0.333), Some(0.667), None]);
        t(&mut buf, &array, b"0.333,0.667,.")?;

        let array = ArrayBuf::Character(vec![Some('n')]);
        t(&mut buf, &array, b"n")?;

        let array = ArrayBuf::Character(vec![Some('n'), Some(';'), Some('='), Some(':'), None]);
        t(&mut buf, &array, b"n,;,=,%3A,.")?;

        let array = ArrayBuf::String(vec![Some(String::from("noodles"))]);
        t(&mut buf, &array, b"noodles")?;

        let array = ArrayBuf::String(vec![
            Some(String::from("noodles")),
            Some(String::from(";")),
            Some(String::from("=")),
            Some(String::from(":")),
            None,
        ]);
        t(&mut buf, &array, b"noodles,;,=,%3A,.")?;

        Ok(())
    }
}
