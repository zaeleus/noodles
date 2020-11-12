#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FileDefinition {
    major_version: u8,
    minor_version: u8,
    file_id: [u8; 20],
}

impl FileDefinition {
    pub fn new(major_version: u8, minor_version: u8, file_id: [u8; 20]) -> Self {
        Self {
            major_version,
            minor_version,
            file_id,
        }
    }

    pub fn major_version(&self) -> u8 {
        self.major_version
    }

    pub fn minor_version(&self) -> u8 {
        self.minor_version
    }

    pub fn file_id(&self) -> &[u8] {
        &self.file_id
    }
}
