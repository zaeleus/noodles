use std::{
    error, fmt,
    io::{self, Write},
};

use super::MISSING;
use crate::variant::record::AlternateBases;

/// An error returns when alternate bases fail to write.
#[derive(Debug)]
pub enum WriteError {
    /// I/O error.
    Io(io::Error),
    /// An alternate base(s) is invalid.
    InvalidAlternateBases(String),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidAlternateBases(_) => None,
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidAlternateBases(s) => write!(f, "invalid alternate base(s): {s}"),
        }
    }
}

pub(super) fn write_alternate_bases<W, B>(
    writer: &mut W,
    alternate_bases: B,
) -> Result<(), WriteError>
where
    W: Write,
    B: AlternateBases,
{
    if alternate_bases.is_empty() {
        writer.write_all(MISSING).map_err(WriteError::Io)?;
    } else {
        for (i, result) in alternate_bases.iter().enumerate() {
            let bases = result.map_err(WriteError::Io)?;

            if i > 0 {
                write!(writer, ",").map_err(WriteError::Io)?;
            }

            if is_valid(bases) {
                write!(writer, "{bases}").map_err(WriteError::Io)?;
            } else {
                return Err(WriteError::InvalidAlternateBases(bases.into()));
            }
        }
    }

    Ok(())
}

// § 1.6.1.5 "Fixed fields: ALT" (2023-08-23): "...no whitespace, commas..."
fn is_valid(s: &str) -> bool {
    fn is_valid_char(c: char) -> bool {
        !c.is_whitespace() && c != ','
    }

    s.chars().all(is_valid_char)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record_buf::AlternateBases as AlternateBasesBuf;

    #[test]
    fn test_write_alternate_bases() -> Result<(), WriteError> {
        let mut buf = Vec::new();

        buf.clear();
        let alternate_bases = AlternateBasesBuf::default();
        write_alternate_bases(&mut buf, &alternate_bases)?;
        assert_eq!(buf, b".");

        buf.clear();
        let alternate_bases = AlternateBasesBuf::from(vec![String::from("C")]);
        write_alternate_bases(&mut buf, &alternate_bases)?;
        assert_eq!(buf, b"C");

        buf.clear();
        let alternate_bases = AlternateBasesBuf::from(vec![String::from("C"), String::from("GT")]);
        write_alternate_bases(&mut buf, &alternate_bases)?;
        assert_eq!(buf, b"C,GT");

        buf.clear();
        let alternate_bases = AlternateBasesBuf::from(vec![String::from("C,")]);
        assert!(matches!(
            write_alternate_bases(&mut buf, &alternate_bases),
            Err(WriteError::InvalidAlternateBases(s)) if s == "C,"
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid("C"));
        assert!(is_valid("GT"));
        assert!(is_valid("*"));
        assert!(is_valid("."));
        assert!(is_valid("<DEL>"));
        assert!(is_valid("<*>"));
        assert!(is_valid("<NON_REF>"));
        assert!(is_valid("C[1:1["));
        assert!(is_valid("]1:0]A"));

        assert!(!is_valid("C,"));
        assert!(!is_valid("G T"));
    }
}
