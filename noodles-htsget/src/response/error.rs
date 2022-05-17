use std::{error, fmt};

use serde::Deserialize;

/// The error kind.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq)]
pub enum Kind {
    // The authorization provided is invalid.
    InvalidAuthentication,
    /// Authorization is required to access the resource.
    PermissionDenied,
    /// The resource was not found.
    NotFound,
    /// The request size is too large.
    PayloadTooLarge,
    /// The requested file format is not supported by the server.
    UnsupportedFormat,
    /// The request parameters are invalid.
    InvalidInput,
    /// The request range is invalid.
    InvalidRange,
}

/// An htsget error.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
pub struct Error {
    #[serde(rename = "error")]
    kind: Kind,
    message: String,
}

impl error::Error for Error {}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}: {}", self.kind, self.message)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let error = Error {
            kind: Kind::NotFound,
            message: String::from("The requested resource was not found"),
        };

        assert_eq!(
            error.to_string(),
            "NotFound: The requested resource was not found"
        );
    }
}
