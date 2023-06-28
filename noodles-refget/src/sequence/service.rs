mod builder;

pub use self::builder::Builder;

use serde::Deserialize;

/// Service information.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
pub struct Service {
    circular_supported: bool,
    algorithms: Vec<String>,
    subsequence_limit: Option<u32>,
    supported_api_version: Vec<String>,
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

    /// Returns the maximum length of an interval.
    ///
    /// If missing, there is no limit.
    pub fn subsequence_limit(&self) -> Option<u32> {
        self.subsequence_limit
    }

    /// Returns a list of supported refget versions.
    pub fn supported_api_versions(&self) -> &[String] {
        &self.supported_api_version
    }
}
