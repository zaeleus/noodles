//! VCF header file format.

/// A VCF header file format.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct FileFormat {
    major: u32,
    minor: u32,
}

const MAJOR_VERSION: u32 = 4;
const MINOR_VERSION: u32 = 4;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let file_format = FileFormat::default();
        assert_eq!(file_format.major(), MAJOR_VERSION);
        assert_eq!(file_format.minor(), MINOR_VERSION);
    }
}
