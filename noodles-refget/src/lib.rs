#![warn(missing_docs)]

//! **noodles-refget** is a refget client.

mod client;
mod sequence;

pub use self::{client::Client, sequence::Sequence};

use std::{error, fmt};

type Result<T> = std::result::Result<T, Error>;

/// An error returned when anything fails to process.
#[derive(Debug)]
pub enum Error {
    /// The URL failed to parse.
    Url(url::ParseError),
    /// The request failed to process.
    Request(reqwest::Error),
}

impl error::Error for Error {}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Url(e) => write!(f, "URL error: {}", e),
            Self::Request(e) => write!(f, "request error: {}", e),
        }
    }
}
