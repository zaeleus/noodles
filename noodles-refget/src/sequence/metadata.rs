mod alias;

pub use self::alias::Alias;

use serde::Deserialize;

/// Sequence metadata.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
pub struct Metadata {
    md5: String,
    trunc512: Option<String>,
    length: u32,
    aliases: Vec<Alias>,
}

impl Metadata {
    /// Returns the MD5 digest in hexadecimal.
    pub fn md5(&self) -> &str {
        &self.md5
    }

    /// Returns the TRUNC512 digest in hexadecimal.
    pub fn trunc512(&self) -> Option<&str> {
        self.trunc512.as_deref()
    }

    /// Returns the length.
    pub fn length(&self) -> u32 {
        self.length
    }

    /// Returns the known aliases.
    pub fn aliases(&self) -> &[Alias] {
        &self.aliases
    }
}
