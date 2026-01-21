use std::{
    error, fmt,
    io::{self, Write},
};

/// An error returns when an info field key fails to write.
#[derive(Debug)]
pub enum WriteError {
    // I/O error.
    Io(io::Error),
    /// The input is invalid.
    Invalid(String),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::Invalid(_) => None,
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::Invalid(s) => write!(f, "invalid input: {s}"),
        }
    }
}

pub(super) fn write_key<W>(writer: &mut W, key: &str) -> Result<(), WriteError>
where
    W: Write,
{
    if is_valid(key) {
        writer.write_all(key.as_bytes()).map_err(WriteError::Io)
    } else {
        Err(WriteError::Invalid(key.into()))
    }
}

// § 1.6.1.8 "Fixed fields: INFO" (2023-08-23): "`^([A-Za-z_][0-9A-Za-z_.]*|1000G)$`".
fn is_valid(s: &str) -> bool {
    use crate::variant::record::info::field::key;

    fn is_valid_char(c: char) -> bool {
        c.is_ascii_alphanumeric() || matches!(c, '_' | '.')
    }

    let mut chars = s.chars();

    let is_valid_first_char = chars
        .next()
        .map(|c| c.is_ascii_alphabetic() || c == '_')
        .unwrap_or_default();

    (is_valid_first_char && chars.all(is_valid_char)) || s == key::IS_IN_1000_GENOMES
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_key() -> Result<(), WriteError> {
        let mut buf = Vec::new();

        buf.clear();
        write_key(&mut buf, "NS")?;
        assert_eq!(buf, b"NS");

        assert!(matches!(
            write_key(&mut buf, "A A"),
            Err(WriteError::Invalid(s)) if s == "A A"
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid("NS"));
        assert!(is_valid("MQ0"));
        assert!(is_valid("1000G"));

        assert!(!is_valid(""));
        assert!(!is_valid("A A"));
    }
}
