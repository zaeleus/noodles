use std::{
    error, fmt,
    io::{self, Write},
};

/// An error returns when a sample integer value fails to write.
#[derive(Debug)]
pub enum WriteError {
    // I/O error.
    Io(io::Error),
    /// The input is invalid.
    Invalid(i32),
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
            Self::Invalid(n) => write!(f, "invalid input: expected (-2^31 + 7) < {n} < 2^31"),
        }
    }
}

pub(super) fn write_integer<W>(writer: &mut W, n: i32) -> Result<(), WriteError>
where
    W: Write,
{
    if is_valid(n) {
        write!(writer, "{n}").map_err(WriteError::Io)
    } else {
        Err(WriteError::Invalid(n))
    }
}

// § 1.3 "Data types" (2024-06-28): "For the Integer type, the values from -2^31 to -2^31 + 7
// cannot be stored in the binary version and therefore are disallowed in both VCF and BCF..."
fn is_valid(n: i32) -> bool {
    const MIN: i32 = i32::MIN + 7;
    n > MIN
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_integer() -> Result<(), WriteError> {
        fn t(buf: &mut Vec<u8>, n: i32, expected: &[u8]) -> Result<(), WriteError> {
            buf.clear();
            write_integer(buf, n)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, i32::MIN + 8, b"-2147483640")?;
        t(&mut buf, 0, b"0")?;
        t(&mut buf, i32::MAX, b"2147483647")?;

        buf.clear();
        assert!(matches!(
            write_integer(&mut buf, i32::MIN),
            Err(WriteError::Invalid(i32::MIN))
        ));

        Ok(())
    }
}
