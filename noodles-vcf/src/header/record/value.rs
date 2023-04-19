//! VCF header record value.

mod collection;
pub mod map;

pub use self::{collection::Collection, map::Map};

/// A VCF header record other value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Value {
    /// A string.
    String(String),
    /// A map.
    Map(String, Map<map::Other>),
}

impl From<&str> for Value {
    fn from(s: &str) -> Self {
        Self::String(s.into())
    }
}

impl From<String> for Value {
    fn from(s: String) -> Self {
        Self::String(s)
    }
}

impl From<(String, Map<map::Other>)> for Value {
    fn from((id, map): (String, Map<map::Other>)) -> Self {
        Self::Map(id, map)
    }
}
