/// A CRAM file definition version.
///
/// This is also called the format number.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Version {
    major: u8,
    minor: u8,
}

impl Version {
    /// Creates a file definition version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::file_definition::Version;
    /// let version = Version::new(3, 0);
    /// ```
    pub fn new(major: u8, minor: u8) -> Self {
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
}

impl Default for Version {
    fn default() -> Self {
        Self::new(3, 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(Version::default(), Version::new(3, 0));
    }
}
