//! VCF header file format.

use std::{error, fmt, num, str::FromStr};

/// A VCF header file format.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct FileFormat {
    major: u32,
    minor: u32,
}

static PREFIX: &str = "VCFv";

const MAJOR_VERSION: u32 = 4;
const MINOR_VERSION: u32 = 4;

const DELIMITER: char = '.';

impl FileFormat {
    /// Creates a file format.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::FileFormat;
    /// let file_format = FileFormat::new(4, 3);
    /// ```
    pub const fn new(major: u32, minor: u32) -> Self {
        Self { major, minor }
    }

    /// Returns the major version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::FileFormat;
    /// let file_format = FileFormat::new(4, 3);
    /// assert_eq!(file_format.major(), 4);
    /// ```
    pub const fn major(&self) -> u32 {
        self.major
    }

    /// Returns the minor version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::FileFormat;
    /// let file_format = FileFormat::new(4, 3);
    /// assert_eq!(file_format.minor(), 3);
    /// ```
    pub const fn minor(&self) -> u32 {
        self.minor
    }
}

impl Default for FileFormat {
    fn default() -> Self {
        Self {
            major: MAJOR_VERSION,
            minor: MINOR_VERSION,
        }
    }
}

impl fmt::Display for FileFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}{}{}", PREFIX, self.major(), DELIMITER, self.minor())
    }
}

/// An error returned when a raw VCF header file format fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The prefix is invalid.
    InvalidPrefix,
    /// The version is missing.
    MissingVersion,
    /// The version is invalid.
    InvalidVersion,
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
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidPrefix => f.write_str("invalid prefix"),
            Self::MissingVersion => f.write_str("missing version"),
            Self::InvalidVersion => f.write_str("invalid version"),
            Self::InvalidMajorVersion(_) => f.write_str("invalid major version"),
            Self::InvalidMinorVersion(_) => f.write_str("invalid minor version"),
        }
    }
}

impl FromStr for FileFormat {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let raw_version = s.strip_prefix(PREFIX).ok_or(ParseError::InvalidPrefix)?;

        if raw_version.is_empty() {
            return Err(ParseError::MissingVersion);
        }

        let (raw_major, raw_minor) = raw_version
            .split_once(DELIMITER)
            .ok_or(ParseError::InvalidVersion)?;

        let major = raw_major.parse().map_err(ParseError::InvalidMajorVersion)?;
        let minor = raw_minor.parse().map_err(ParseError::InvalidMinorVersion)?;

        Ok(Self::new(major, minor))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let file_format = FileFormat::default();
        assert_eq!(file_format.major(), MAJOR_VERSION);
        assert_eq!(file_format.minor(), MINOR_VERSION);
    }

    #[test]
    fn test_fmt() {
        let file_format = FileFormat::new(4, 3);
        assert_eq!(file_format.to_string(), "VCFv4.3");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("VCFv4.3".parse(), Ok(FileFormat::new(4, 3)));

        assert_eq!("".parse::<FileFormat>(), Err(ParseError::Empty));

        assert_eq!("4.3".parse::<FileFormat>(), Err(ParseError::InvalidPrefix));
        assert_eq!(
            "NDLv4.3".parse::<FileFormat>(),
            Err(ParseError::InvalidPrefix)
        );

        assert_eq!(
            "VCFv".parse::<FileFormat>(),
            Err(ParseError::MissingVersion)
        );

        assert_eq!(
            "VCFvx".parse::<FileFormat>(),
            Err(ParseError::InvalidVersion)
        );

        assert_eq!(
            "VCFv4".parse::<FileFormat>(),
            Err(ParseError::InvalidVersion)
        );

        assert!(matches!(
            "VCFvx.3".parse::<FileFormat>(),
            Err(ParseError::InvalidMajorVersion(_))
        ));

        assert!(matches!(
            "VCFv4.x".parse::<FileFormat>(),
            Err(ParseError::InvalidMinorVersion(_))
        ));
    }
}
