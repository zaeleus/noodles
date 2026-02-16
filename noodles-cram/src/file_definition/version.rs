use std::{cmp::Ordering, fmt, io};

/// A CRAM file definition version.
///
/// This is also called the format number.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Version {
    major: u8,
    minor: u8,
}

impl Version {
    /// CRAM 2.0
    pub const V2_0: Self = Self::new(2, 0);

    /// CRAM 2.1
    pub const V2_1: Self = Self::new(2, 1);

    /// CRAM 3.0
    pub const V3_0: Self = Self::new(3, 0);

    /// CRAM 3.1
    pub const V3_1: Self = Self::new(3, 1);

    /// CRAM 4.0
    pub const V4_0: Self = Self::new(4, 0);

    /// Creates a file definition version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::file_definition::Version;
    /// let version = Version::new(3, 0);
    /// ```
    pub const fn new(major: u8, minor: u8) -> Self {
        Self { major, minor }
    }

    /// Returns the major version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::file_definition::Version;
    /// let version = Version::new(3, 0);
    /// assert_eq!(version.major(), 3);
    /// ```
    pub fn major(&self) -> u8 {
        self.major
    }

    /// Returns the minor version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::file_definition::Version;
    /// let version = Version::new(3, 0);
    /// assert_eq!(version.minor(), 0);
    /// ```
    pub fn minor(&self) -> u8 {
        self.minor
    }

    /// Returns `true` if this version uses CRC32 checksums (>= 3.0).
    pub fn has_crc32(&self) -> bool {
        *self >= Self::V3_0
    }

    /// Returns `true` if this version uses VLQ (uint7/sint7) varints instead of ITF8/LTF8 (>= 4.0).
    pub fn uses_vlq(&self) -> bool {
        *self >= Self::V4_0
    }

    /// Returns `true` if this version uses 64-bit positions (>= 4.0).
    pub fn has_64bit_positions(&self) -> bool {
        *self >= Self::V4_0
    }

    /// Validates that the version is a supported CRAM version for reading.
    pub fn validate(&self) -> io::Result<()> {
        match *self {
            Self::V2_0 | Self::V2_1 | Self::V3_0 | Self::V3_1 | Self::V4_0 => Ok(()),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("unsupported CRAM version: {self}"),
            )),
        }
    }
}

impl Default for Version {
    fn default() -> Self {
        Self::new(3, 0)
    }
}

impl PartialOrd for Version {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Version {
    fn cmp(&self, other: &Self) -> Ordering {
        self.major
            .cmp(&other.major)
            .then(self.minor.cmp(&other.minor))
    }
}

impl fmt::Display for Version {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}.{}", self.major, self.minor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(Version::default(), Version::new(3, 0));
    }

    #[test]
    fn test_ordering() {
        assert!(Version::V2_0 < Version::V2_1);
        assert!(Version::V2_1 < Version::V3_0);
        assert!(Version::V3_0 < Version::V3_1);
        assert!(Version::V3_1 < Version::V4_0);
    }

    #[test]
    fn test_has_crc32() {
        assert!(!Version::V2_0.has_crc32());
        assert!(!Version::V2_1.has_crc32());
        assert!(Version::V3_0.has_crc32());
        assert!(Version::V3_1.has_crc32());
        assert!(Version::V4_0.has_crc32());
    }

    #[test]
    fn test_uses_vlq() {
        assert!(!Version::V2_0.uses_vlq());
        assert!(!Version::V3_0.uses_vlq());
        assert!(!Version::V3_1.uses_vlq());
        assert!(Version::V4_0.uses_vlq());
    }

    #[test]
    fn test_has_64bit_positions() {
        assert!(!Version::V3_0.has_64bit_positions());
        assert!(Version::V4_0.has_64bit_positions());
    }

    #[test]
    fn test_validate() {
        assert!(Version::V2_0.validate().is_ok());
        assert!(Version::V2_1.validate().is_ok());
        assert!(Version::V3_0.validate().is_ok());
        assert!(Version::V3_1.validate().is_ok());
        assert!(Version::V4_0.validate().is_ok());
        assert!(Version::new(1, 0).validate().is_err());
        assert!(Version::new(5, 0).validate().is_err());
    }
}
