//! VCF header record other key.

use std::{borrow::Borrow, fmt};

/// A nonstandard VCF record key.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Other(pub(super) String);

impl AsRef<str> for Other {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl Borrow<str> for Other {
    fn borrow(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for Other {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.as_ref().fmt(f)
    }
}
