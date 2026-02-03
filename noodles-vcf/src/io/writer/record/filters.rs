use std::{
    error, fmt,
    io::{self, Write},
};

use super::MISSING;
use crate::{Header, variant::record::Filters};

/// An error returns when filters fail to write.
#[derive(Debug)]
pub enum WriteError {
    /// I/O error.
    Io(io::Error),
    /// A filter is invalid.
    InvalidFilter(String),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidFilter(_) => None,
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidFilter(s) => write!(f, "invalid filter: {s}"),
        }
    }
}

pub(super) fn write_filters<W, F>(
    writer: &mut W,
    header: &Header,
    filters: F,
) -> Result<(), WriteError>
where
    W: Write,
    F: Filters,
{
    const DELIMITER: &[u8] = b";";

    if filters.is_empty() {
        writer.write_all(MISSING).map_err(WriteError::Io)?;
    } else {
        for (i, result) in filters.iter(header).enumerate() {
            let id = result.map_err(WriteError::Io)?;

            if i > 0 {
                writer.write_all(DELIMITER).map_err(WriteError::Io)?;
            }

            if is_valid(id) {
                writer.write_all(id.as_bytes()).map_err(WriteError::Io)?;
            } else {
                return Err(WriteError::InvalidFilter(id.into()));
            }
        }
    }

    Ok(())
}

// § 1.6.1.7 "Fixed fields: FILTER" (2023-08-23): "...no whitespace or semicolons permitted..."
fn is_valid(s: &str) -> bool {
    fn is_valid_char(c: char) -> bool {
        !c.is_whitespace() && c != ';'
    }

    s.chars().all(is_valid_char)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record_buf::Filters as FiltersBuf;

    #[test]
    fn test_write_filters() -> Result<(), WriteError> {
        fn t(
            buf: &mut Vec<u8>,
            header: &Header,
            filters: &FiltersBuf,
            expected: &[u8],
        ) -> Result<(), WriteError> {
            buf.clear();
            write_filters(buf, header, filters)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let header = Header::default();
        let mut buf = Vec::new();

        t(&mut buf, &header, &FiltersBuf::default(), b".")?;
        t(&mut buf, &header, &FiltersBuf::pass(), b"PASS")?;

        let filters = [String::from("q10")].into_iter().collect();
        t(&mut buf, &header, &filters, b"q10")?;

        let filters = [String::from("q10"), String::from("s50")]
            .into_iter()
            .collect();
        t(&mut buf, &header, &filters, b"q10;s50")?;

        buf.clear();
        let filters = [String::from("q 10")].into_iter().collect();
        assert!(matches!(
            write_filters(&mut buf, &header, &filters),
            Err(WriteError::InvalidFilter(s)) if s == "q 10"
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid("PASS"));
        assert!(is_valid("q10"));

        assert!(!is_valid("q 10"));
        assert!(!is_valid("q10;"));
    }
}
