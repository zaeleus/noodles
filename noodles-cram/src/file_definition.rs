mod version;

pub use self::version::Version;

#[derive(Clone, Default, Debug, Eq, PartialEq)]
pub struct FileDefinition {
    version: Version,
    file_id: [u8; 20],
}

impl FileDefinition {
    pub fn new(version: Version, file_id: [u8; 20]) -> Self {
        Self { version, file_id }
    }

    pub fn version(&self) -> Version {
        self.version
    }

    pub fn file_id(&self) -> &[u8] {
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
