mod builder;

pub use self::builder::Builder;

use serde::Deserialize;

/// Service information.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
pub struct Service {
    circular_supported: bool,
    algorithms: Vec<String>,
    identifier_types: Vec<String>,
    subsequence_limit: Option<u32>,
}

impl Service {
    /// Returns whether circular genomes are supported by the server.
    pub fn circular_supported(&self) -> bool {
        self.circular_supported
    }

    /// Returns a list of digest algorithms supported by the server.
    pub fn algorithms(&self) -> &[String] {
        &self.algorithms
    }

    /// Returns a list of supported sequence type identifiers.
    pub fn identifier_types(&self) -> &[String] {
        &self.identifier_types
    }

    /// Returns the maximum length of an interval.
    ///
    /// If missing or smaller than 1, there is no limit.
    pub fn subsequence_limit(&self) -> Option<u32> {
        self.subsequence_limit
    }
}
