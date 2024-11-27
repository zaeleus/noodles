//! GFF directive version.

use std::{error, fmt, num, str::FromStr};

const MAJOR_VERSION: u32 = 3;

const DELIMITER: char = '.';
const MAX_COMPONENT_COUNT: usize = 3;

/// A GFF directive version.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GffVersion {
    major: u32,
    minor: Option<u32>,
    patch: Option<u32>,
}

impl GffVersion {
    /// Returns the major version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::directive_buf::value::GffVersion;
    /// let version: GffVersion = "3.1.26".parse()?;
    /// assert_eq!(version.major(), 3);
    /// # Ok::<(), noodles_gff::directive_buf::value::gff_version::ParseError>(())
    /// ```
    pub fn major(&self) -> u32 {
        self.major
    }

    /// Returns the minor version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::directive_buf::value::GffVersion;
    /// let version: GffVersion = "3.1.26".parse()?;
    /// assert_eq!(version.minor(), Some(1));
    /// # Ok::<(), noodles_gff::directive_buf::value::gff_version::ParseError>(())
    /// ```
    pub fn minor(&self) -> Option<u32> {
        self.minor
    }

    /// Returns the patch version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::directive_buf::value::GffVersion;
    /// let version: GffVersion = "3.1.26".parse()?;
    /// assert_eq!(version.patch(), Some(26));
    /// # Ok::<(), noodles_gff::directive_buf::value::gff_version::ParseError>(())
    /// ```
    pub fn patch(&self) -> Option<u32> {
        self.patch
    }
}

impl Default for GffVersion {
    fn default() -> Self {
        Self {
            major: MAJOR_VERSION,
            minor: None,
            patch: None,
        }
    }
}

impl fmt::Display for GffVersion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.major())?;

        if let Some(minor) = self.minor() {
            write!(f, "{DELIMITER}{minor}")?;

            if let Some(patch) = self.patch() {
                write!(f, "{DELIMITER}{patch}")?;
            }
        }

        Ok(())
    }
}

/// An error returned when a raw GFF version directive fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The major version is missing.
    MissingMajorVersion,
    /// The major version is invalid.
    InvalidMajorVersion(num::ParseIntError),
    /// The minor version is invalid.
    InvalidMinorVersion(num::ParseIntError),
    /// The patch version is invalid.
    InvalidPatchVersion(num::ParseIntError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidMajorVersion(e)
            | Self::InvalidMinorVersion(e)
            | Self::InvalidPatchVersion(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::MissingMajorVersion => f.write_str("missing major version"),
            Self::InvalidMajorVersion(_) => f.write_str("invalid major version"),
            Self::InvalidMinorVersion(_) => f.write_str("invalid minor version"),
            Self::InvalidPatchVersion(_) => f.write_str("invalid patch version"),
        }
    }
}

impl FromStr for GffVersion {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let mut components = s.splitn(MAX_COMPONENT_COUNT, DELIMITER);

        let major = components
            .next()
            .ok_or(ParseError::MissingMajorVersion)
            .and_then(|t| t.parse().map_err(ParseError::InvalidMajorVersion))?;

        let minor = match components.next() {
            Some(t) => t
                .parse()
                .map(Some)
                .map_err(ParseError::InvalidMinorVersion)?,
            None => None,
        };

        let patch = match components.next() {
            Some(t) => t
                .parse()
                .map(Some)
                .map_err(ParseError::InvalidPatchVersion)?,
            None => None,
        };

        Ok(Self {
            major,
            minor,
            patch,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let version = GffVersion {
            major: 3,
            minor: None,
            patch: None,
        };
        assert_eq!(version.to_string(), "3");

        let version = GffVersion {
            major: 3,
            minor: Some(1),
            patch: None,
        };
        assert_eq!(version.to_string(), "3.1");

        let version = GffVersion {
            major: 3,
            minor: Some(1),
            patch: Some(26),
        };
        assert_eq!(version.to_string(), "3.1.26");
    }

    #[test]
    fn test_from_str() {
        assert_eq!(
            "3".parse(),
            Ok(GffVersion {
                major: 3,
                minor: None,
                patch: None
            })
        );

        assert_eq!(
            "3.1".parse(),
            Ok(GffVersion {
                major: 3,
                minor: Some(1),
                patch: None
            })
        );

        assert_eq!(
            "3.1.26".parse(),
            Ok(GffVersion {
                major: 3,
                minor: Some(1),
                patch: Some(26),
            })
        );

        assert_eq!("".parse::<GffVersion>(), Err(ParseError::Empty));

        assert!(matches!(
            "a".parse::<GffVersion>(),
            Err(ParseError::InvalidMajorVersion(_))
        ));

        assert!(matches!(
            "3.b".parse::<GffVersion>(),
            Err(ParseError::InvalidMinorVersion(_))
        ));

        assert!(matches!(
            "3.1.c".parse::<GffVersion>(),
            Err(ParseError::InvalidPatchVersion(_))
        ));
    }
}
