mod builder;

pub use self::builder::Builder;

use std::num::NonZeroU32;

use serde::{Deserialize, Deserializer};

/// Service information.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
pub struct Service {
    circular_supported: bool,
    algorithms: Vec<String>,
    identifier_types: Vec<String>,
    #[serde(deserialize_with = "deserialize_subsequence_limit")]
    subsequence_limit: Option<NonZeroU32>,
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
    /// If missing, there is no limit.
    pub fn subsequence_limit(&self) -> Option<NonZeroU32> {
        self.subsequence_limit
    }
}

fn deserialize_subsequence_limit<'de, D>(deserializer: D) -> Result<Option<NonZeroU32>, D::Error>
where
    D: Deserializer<'de>,
{
    let value: Option<u32> = Deserialize::deserialize(deserializer)?;

    match value {
        Some(n) => Ok(NonZeroU32::new(n)),
        None => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use serde_test::{Token, assert_de_tokens};

    use super::*;

    #[test]
    fn test_deserialize_subsequence_limit() {
        #[derive(Debug, Deserialize, Eq, PartialEq)]
        #[serde(transparent)]
        struct SubsequenceLimit(
            #[serde(deserialize_with = "deserialize_subsequence_limit")] Option<NonZeroU32>,
        );

        assert_de_tokens(
            &SubsequenceLimit(NonZeroU32::new(8)),
            &[Token::Some, Token::I64(8)],
        );
        assert_de_tokens(&SubsequenceLimit(None), &[Token::Some, Token::I64(0)]);
        assert_de_tokens(&SubsequenceLimit(None), &[Token::None]);
    }
}
