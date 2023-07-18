mod alias;
pub(crate) mod builder;

pub use self::alias::Alias;
pub(crate) use self::builder::Builder;

use serde::Deserialize;

/// Sequence metadata.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
pub struct Metadata {
    md5: String,
    ga4gh: Option<String>,
    length: u32,
    aliases: Vec<Alias>,
}

impl Metadata {
    /// Returns the MD5 digest in hexadecimal.
    pub fn md5(&self) -> &str {
        &self.md5
    }

    /// Returns the ga4gh digest in hexadecimal.
    pub fn ga4gh(&self) -> Option<&str> {
        self.ga4gh.as_deref()
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
