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
            // SAM 1.6 no longer allows duplicate tags. See § 1.3 "The header section" (2025-08-12):
            // "Within each (non-`@CO`) header line, no field tag may appear more than once..."
            allow_duplicate_tags: version < Version::new(1, 6),
        }
    }
}
