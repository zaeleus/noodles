use std::{
    error, fmt,
    io::{self, Write},
};

use crate::variant::record::samples::keys::key;

/// An error returns when record samples keys fail to write.
#[derive(Debug)]
pub enum WriteError {
    // I/O error.
    Io(io::Error),
    // The genotype (GT) key is not first.
    InvalidGenotypePosition(usize),
    // A key is invalid.
    InvalidKey(String),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidGenotypePosition(_) => None,
            Self::InvalidKey(_) => None,
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidGenotypePosition(i) => {
                write!(f, "invalid genotype (GT) position: expected 0, got {i}")
            }
            Self::InvalidKey(s) => write!(f, "invalid key: {s}"),
        }
    }
}

pub(super) fn write_keys<'a, W, I>(writer: &mut W, keys: I) -> Result<(), WriteError>
where
    W: Write,
    I: Iterator<Item = io::Result<&'a str>>,
{
    const DELIMITER: &[u8] = b":";

    for (i, result) in keys.enumerate() {
        let key = result.map_err(WriteError::Io)?;

        if i > 0 {
            if key == key::GENOTYPE {
                return Err(WriteError::InvalidGenotypePosition(i));
            }

            writer.write_all(DELIMITER).map_err(WriteError::Io)?;
        }

        write_key(writer, key)?;
    }

    Ok(())
}

fn write_key<W>(writer: &mut W, key: &str) -> Result<(), WriteError>
where
    W: Write,
{
    if is_valid(key) {
        writer.write_all(key.as_bytes()).map_err(WriteError::Io)
    } else {
        Err(WriteError::InvalidKey(key.into()))
    }
}

// § 1.6.2 "Genotype fields" (2023-08-23): "`^[A-Za-z_][0-9A-Za-z_.]*$`".
fn is_valid(s: &str) -> bool {
    fn is_valid_char(c: char) -> bool {
        c.is_ascii_alphanumeric() || matches!(c, '_' | '.')
    }

    let mut chars = s.chars();

    let is_valid_first_char = chars
        .next()
        .map(|c| c.is_ascii_alphabetic() || c == '_')
        .unwrap_or_default();

    is_valid_first_char && chars.all(is_valid_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_keys() -> Result<(), WriteError> {
        let mut buf = Vec::new();

        buf.clear();
        let keys = [Ok("GT")];
        write_keys(&mut buf, keys.into_iter())?;
        assert_eq!(buf, b"GT");

        buf.clear();
        let keys = [Ok("GT"), Ok("GQ")];
        write_keys(&mut buf, keys.into_iter())?;
        assert_eq!(buf, b"GT:GQ");

        buf.clear();
        let keys = [Ok("GQ")];
        write_keys(&mut buf, keys.into_iter())?;
        assert_eq!(buf, b"GQ");

        buf.clear();
        let keys = [Ok("GQ"), Ok("GT")];
        assert!(matches!(
            write_keys(&mut buf, keys.into_iter()),
            Err(WriteError::InvalidGenotypePosition(1))
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid("GT"));
        assert!(is_valid("PSL"));

        assert!(!is_valid(""));
        assert!(!is_valid("G T"));
        assert!(!is_valid("1000G"));
    }
}
