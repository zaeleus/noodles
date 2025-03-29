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
    /// The response had an unsuccessful HTTP status code.  
    Response(reqwest::Error),
}

impl error::Error for Error {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Input => None,
            Self::Url(e) => Some(e),
            Self::Request(e) => Some(e),
            Self::Response(e) => Some(e),
        }
    }
}

impl From<reqwest::Error> for Error {
    fn from(err: reqwest::Error) -> Self {
        if err.is_status() {
            Error::Response(err)
        } else {
            Error::Request(err)
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Input => f.write_str("invalid input"),
            Self::Url(_) => f.write_str("URL error"),
            Self::Request(_) => f.write_str("request error"),
            Self::Response(_) => f.write_str("response error"),
        }
    }
}
