use crate::header::record::value::map::header::Version;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Context {
    allow_duplicate_tags: bool,
}

impl Context {
    pub fn allow_duplicate_tags(&self) -> bool {
        self.allow_duplicate_tags
    }
}

impl Default for Context {
    fn default() -> Self {
        Self::from(Version::default())
    }
}

impl From<Version> for Context {
    fn from(version: Version) -> Self {
        Self {
            allow_duplicate_tags: version < Version::new(1, 6),
        }
    }
}
