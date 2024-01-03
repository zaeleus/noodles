//! SAM record name.

use std::{
    error, fmt,
    str::{self, FromStr},
};

/// § 1.4 "The alignment section: mandatory fields" (2021-06-03): "A `QNAME` '*' indicates the
/// information is unavailable".
pub const MISSING: &[u8] = b"*";

// § 1.4 The alignment section: mandatory fields (2020-07-19)
const MAX_LENGTH: usize = 254;

/// A SAM record name.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct Name(Vec<u8>);

impl Name {
    /// Creates a name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::Name;
    /// let name = Name::try_new("r1")?;
    /// # Ok::<_, noodles_sam::record::name::ParseError>(())
    /// ```
    pub fn try_new<I>(data: I) -> Result<Self, ParseError>
    where
        I: Into<Vec<u8>>,
    {
        Self::try_from(data.into())
    }

    /// Returns the length of the name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::Name;
    /// let name: Name = "r1".parse()?;
    /// assert_eq!(name.len(), 2);
    /// # Ok::<_, noodles_sam::record::name::ParseError>(())
    /// ```
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.0.len()
    }
}

impl crate::alignment::record::Name for Name {
    fn as_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl crate::alignment::record::Name for &Name {
    fn as_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl AsRef<[u8]> for Name {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsRef<str> for Name {
    fn as_ref(&self) -> &str {
        // SAFETY: The internal buffer is limited to ASCII graphic characters (i.e., '!'-'~').
        str::from_utf8(&self.0).unwrap()
    }
}

impl fmt::Display for Name {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

/// An error returned when a raw SAM record name fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid => f.write_str("invalid input"),
        }
    }
}

impl FromStr for Name {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::try_new(s)
    }
}

impl TryFrom<Vec<u8>> for Name {
    type Error = ParseError;

    fn try_from(buf: Vec<u8>) -> Result<Self, Self::Error> {
        if buf.is_empty() {
            Err(ParseError::Empty)
        } else if !is_valid_name(&buf) {
            Err(ParseError::Invalid)
        } else {
            Ok(Self(buf))
        }
    }
}

impl From<Name> for Vec<u8> {
    fn from(name: Name) -> Self {
        name.0
    }
}

fn is_valid_name_char(b: u8) -> bool {
    b.is_ascii_graphic() && b != b'@'
}

// § 1.4 "The alignment section: mandatory fields" (2021-06-03): "`[!-?A-~]{1,254}`".
fn is_valid_name(buf: &[u8]) -> bool {
    buf != MISSING && buf.len() <= MAX_LENGTH && buf.iter().copied().all(is_valid_name_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let name: Name = "r0".parse()?;
        assert_eq!(name.to_string(), "r0");
        Ok(())
    }

    #[test]
    fn test_try_new() {
        assert_eq!(Name::try_new("r0"), Ok(Name(b"r0".to_vec())));

        assert_eq!(Name::try_new(""), Err(ParseError::Empty));
        assert_eq!(Name::try_new("*"), Err(ParseError::Invalid));
        assert_eq!(Name::try_new("r 0"), Err(ParseError::Invalid));
        assert_eq!(Name::try_new("@r0"), Err(ParseError::Invalid));

        let s = "n".repeat(MAX_LENGTH + 1);
        assert_eq!(Name::try_new(s), Err(ParseError::Invalid));
    }
}
