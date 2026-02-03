use std::{
    error, fmt,
    io::{self, Write},
};

use super::MISSING;
use crate::variant::record::Ids;

/// An error returns when IDs fails to write.
#[derive(Debug)]
pub enum WriteError {
    /// I/O error.
    Io(io::Error),
    /// An ID is invalid.
    InvalidId(String),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidId(_) => None,
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidId(s) => write!(f, "invalid ID: {s}"),
        }
    }
}

pub(super) fn write_ids<W, I>(writer: &mut W, ids: I) -> Result<(), WriteError>
where
    W: Write,
    I: Ids,
{
    const DELIMITER: &[u8] = b";";

    if ids.is_empty() {
        writer.write_all(MISSING).map_err(WriteError::Io)?;
    } else {
        for (i, id) in ids.iter().enumerate() {
            if i > 0 {
                writer.write_all(DELIMITER).map_err(WriteError::Io)?;
            }

            if is_valid(id) {
                writer.write_all(id.as_bytes()).map_err(WriteError::Io)?;
            } else {
                return Err(WriteError::InvalidId(id.into()));
            }
        }
    }

    Ok(())
}

// § 1.6.1.3 "Fixed fields: ID" (2023-08-23): "...no whitespace or semicolons permitted..."
fn is_valid(s: &str) -> bool {
    fn is_valid_char(c: char) -> bool {
        !c.is_whitespace() && c != ';'
    }

    s.chars().all(is_valid_char)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record_buf::Ids as IdsBuf;

    #[test]
    fn test_write_ids() -> Result<(), WriteError> {
        fn t(buf: &mut Vec<u8>, ids: &IdsBuf, expected: &[u8]) -> Result<(), WriteError> {
            buf.clear();
            write_ids(buf, ids)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let ids = IdsBuf::default();
        t(&mut buf, &ids, b".")?;

        let ids = [String::from("id0")].into_iter().collect();
        t(&mut buf, &ids, b"id0")?;

        let ids = [String::from("id0"), String::from("id1")]
            .into_iter()
            .collect();
        t(&mut buf, &ids, b"id0;id1")?;

        buf.clear();
        let ids = [String::from("id 0")].into_iter().collect();
        assert!(matches!(
            write_ids(&mut buf, &ids),
            Err(WriteError::InvalidId(s)) if s == "id 0"
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid("id0"));

        assert!(!is_valid("id 0"));
        assert!(!is_valid("id0;"));
    }
}
