//! SAM header header version.

use std::{error, fmt, num, str::FromStr};

const MAJOR_VERSION: u32 = 1;
const MINOR_VERSION: u32 = 6;

const DELIMITER: char = '.';

/// A SAM header header version (`VN`).
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct Version {
    major: u32,
    minor: u32,
}

impl Version {
    /// Creates a new version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::map::header::Version;
    /// let version = Version::new(1, 6);
    /// ```
    pub const fn new(major: u32, minor: u32) -> Self {
        Self { major, minor }
    }

    /// Returns the major version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::map::header::Version;
    /// let version = Version::new(1, 6);
    /// assert_eq!(version.major(), 1);
    /// ```
    pub const fn major(&self) -> u32 {
        self.major
    }

    /// Returns the minor version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::map::header::Version;
    /// let version = Version::new(1, 6);
    /// assert_eq!(version.minor(), 6);
    /// ```
    pub const fn minor(&self) -> u32 {
        self.minor
    }
}

impl Default for Version {
    fn default() -> Self {
        Self::new(MAJOR_VERSION, MINOR_VERSION)
    }
}

impl fmt::Display for Version {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}{}", self.major(), DELIMITER, self.minor())
    }
}

/// An error returned when a raw SAM header header version fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
    /// The major version is invalid.
    InvalidMajorVersion(num::ParseIntError),
    /// The minor version is invalid.
    InvalidMinorVersion(num::ParseIntError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidMajorVersion(e) | Self::InvalidMinorVersion(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    #[allow(deprecated)]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidMajorVersion(_) => f.write_str("invalid major version"),
            Self::InvalidMinorVersion(_) => f.write_str("invalid minor version"),
        }
    }
}

impl FromStr for Version {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        match s.split_once(DELIMITER) {
            Some((a, b)) => {
                let major = a.parse().map_err(ParseError::InvalidMajorVersion)?;
                let minor = b.parse().map_err(ParseError::InvalidMinorVersion)?;
                Ok(Self::new(major, minor))
            }
            None => Err(ParseError::Invalid),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let version = Version::default();
        assert_eq!(version, Version::new(MAJOR_VERSION, MINOR_VERSION));
    }

    #[test]
    fn test_fmt() {
        let version = Version::new(1, 6);
        assert_eq!(version.to_string(), "1.6");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("1.6".parse(), Ok(Version::new(1, 6)));

        assert_eq!("".parse::<Version>(), Err(ParseError::Empty));
        assert_eq!("1".parse::<Version>(), Err(ParseError::Invalid));

        assert!(matches!(
            ".".parse::<Version>(),
            Err(ParseError::InvalidMajorVersion(_))
        ));

        assert!(matches!(
            "x.6".parse::<Version>(),
            Err(ParseError::InvalidMajorVersion(_))
        ));

        assert!(matches!(
            "1.x".parse::<Version>(),
            Err(ParseError::InvalidMinorVersion(_))
        ));

        assert!(matches!(
            "1.6.1".parse::<Version>(),
            Err(ParseError::InvalidMinorVersion(_))
        ));
    }
}
