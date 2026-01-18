use std::{
    error, fmt,
    io::{self, Write},
};

/// An error returns when a reference sequence name fails to write.
#[derive(Debug)]
pub enum WriteError {
    // I/O error.
    Io(io::Error),
    /// The reference sequence name is invalid.
    InvalidReferenceSequenceName(String),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidReferenceSequenceName(_) => None,
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidReferenceSequenceName(s) => {
                write!(f, "invalid reference sequence name: {s}")
            }
        }
    }
}

pub(super) fn write_reference_sequence_name<W>(
    writer: &mut W,
    reference_sequence_name: &str,
) -> Result<(), WriteError>
where
    W: Write,
{
    let name = strip_symbol_delimiters(reference_sequence_name).unwrap_or(reference_sequence_name);

    if is_valid(name) {
        writer
            .write_all(reference_sequence_name.as_bytes())
            .map_err(WriteError::Io)
    } else {
        Err(WriteError::InvalidReferenceSequenceName(
            reference_sequence_name.into(),
        ))
    }
}

fn strip_symbol_delimiters(s: &str) -> Option<&str> {
    const PREFIX: char = '<';
    const SUFFIX: char = '>';

    s.strip_prefix(PREFIX).and_then(|t| t.strip_suffix(SUFFIX))
}

// § 1.4.7 "Contig field format" (2023-08-23): "`[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*`".
fn is_valid(s: &str) -> bool {
    fn is_valid_char(c: char) -> bool {
        ('!'..='~').contains(&c)
            && !matches!(
                c,
                '\\' | ',' | '"' | '`' | '\'' | '(' | ')' | '[' | ']' | '{' | '}' | '<' | '>',
            )
    }

    let mut chars = s.chars();

    let is_valid_first_char = chars
        .next()
        .map(|c| c != '*' && c != '=' && is_valid_char(c))
        .unwrap_or_default();

    is_valid_first_char && chars.all(is_valid_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_reference_sequence_name() -> Result<(), WriteError> {
        fn t(
            buf: &mut Vec<u8>,
            reference_sequence_name: &str,
            expected: &[u8],
        ) -> Result<(), WriteError> {
            buf.clear();
            write_reference_sequence_name(buf, reference_sequence_name)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, "sq0", b"sq0")?;
        t(&mut buf, "<sq0>", b"<sq0>")?;

        buf.clear();
        assert!(matches!(
            write_reference_sequence_name(&mut buf, "sq 0"),
            Err(WriteError::InvalidReferenceSequenceName(s)) if s == "sq 0"
        ));

        buf.clear();
        assert!(matches!(
            write_reference_sequence_name(&mut buf, "<sq 0>"),
            Err(WriteError::InvalidReferenceSequenceName(s)) if s == "<sq 0>"
        ));

        Ok(())
    }

    #[test]
    fn test_strip_symbol_delimiters() {
        assert!(strip_symbol_delimiters("sq0").is_none());
        assert_eq!(strip_symbol_delimiters("<sq0>"), Some("sq0"));
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid("sq0"));

        assert!(!is_valid(""));
        assert!(!is_valid("sq 0"));
        assert!(!is_valid("sq[0]"));
        assert!(!is_valid(">sq0"));
        assert!(!is_valid("*sq0"));
        assert!(!is_valid("=sq0"));
    }
}
