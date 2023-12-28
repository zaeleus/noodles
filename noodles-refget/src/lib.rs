#![warn(missing_docs)]

//! **noodles-refget** is a refget 2.0 client.

mod client;
pub mod sequence;

pub use self::{client::Client, sequence::Sequence};

use std::{error, fmt};

type Result<T> = std::result::Result<T, Error>;

/// An error returned when anything fails to process.
#[derive(Debug)]
pub enum Error {
    /// An input is invalid.
    Input,
    /// The URL failed to parse.
    Url(url::ParseError),
    /// The request failed to process.
    Request(reqwest::Error),
}

impl error::Error for Error {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Input => None,
            Self::Url(e) => Some(e),
            Self::Request(e) => Some(e),
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Input => f.write_str("invalid input"),
            Self::Url(_) => f.write_str("URL error"),
            Self::Request(_) => f.write_str("request error"),
        }
    }
}
