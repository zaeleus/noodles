use std::fmt;

const MAJOR_VERSION: u32 = 1;
const MINOR_VERSION: u32 = 6;

const DELIMITER: char = '.';

/// A SAM header header version (`VN`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
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
    /// use noodles_sam::header::header::Version;
    /// let version = Version::new(1, 6);
    /// ```
    pub fn new(major: u32, minor: u32) -> Self {
        Self { major, minor }
    }

    /// Returns the major version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::Version;
    /// let version = Version::new(1, 6);
    /// assert_eq!(version.major(), 1);
    /// ```
    pub fn major(&self) -> u32 {
        self.major
    }

    /// Returns the minor version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::header::Version;
    /// let version = Version::new(1, 6);
    /// assert_eq!(version.minor(), 6);
    /// ```
    pub fn minor(&self) -> u32 {
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
}
