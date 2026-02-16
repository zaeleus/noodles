//! CRAM file definition and fields.

mod version;

pub use self::version::Version;

/// A CRAM file definition.
///
/// The CRAM file definition holds the format version and file ID. See § 6 File definition
/// (2020-06-22).
#[derive(Clone, Default, Debug, Eq, PartialEq)]
pub struct FileDefinition {
    version: Version,
    file_id: [u8; 20],
}

impl FileDefinition {
    /// Creates a file definition.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::{file_definition::Version, FileDefinition};
    /// let file_definition = FileDefinition::new(Version::new(3, 0), [0; 20]);
    /// ```
    pub fn new(version: Version, file_id: [u8; 20]) -> Self {
        Self { version, file_id }
    }

    /// Returns the file version.
    ///
    /// This is also called the (major and minor) format number.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::{file_definition::Version, FileDefinition};
    /// let file_definition = FileDefinition::new(Version::new(3, 0), [0; 20]);
    /// assert_eq!(file_definition.version(), Version::new(3, 0));
    /// ```
    pub fn version(&self) -> Version {
        self.version
    }

    /// Returns the file ID.
    ///
    /// The file ID has a fixed length of 20 bytes. It can be any arbitrary identifier, e.g., the
    /// file name or a 160-bit checksum.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::{file_definition::Version, FileDefinition};
    /// let file_definition = FileDefinition::new(Version::new(3, 0), [0; 20]);
    /// assert_eq!(file_definition.file_id(), &[0; 20]);
    /// ```
    pub fn file_id(&self) -> &[u8; 20] {
        &self.file_id
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let actual = FileDefinition::default();
        let expected = FileDefinition::new(Version::new(3, 0), [0; 20]);
        assert_eq!(actual, expected);
    }
}
