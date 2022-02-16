use serde::Deserialize;

/// A sequence alias.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
pub struct Alias {
    alias: String,
    naming_authority: String,
}

impl Alias {
    /// Returns the alias.
    pub fn alias(&self) -> &str {
        &self.alias
    }

    /// Returns the naming authority.
    pub fn naming_authority(&self) -> &str {
        &self.naming_authority
    }
}
