#![warn(missing_docs)]

//! **noodles-htsget** is an htsget 1.3 client.

pub(crate) mod chunks;
mod client;
mod format;
pub mod reads;
pub(crate) mod request;
pub mod response;
pub mod variants;

pub use self::{client::Client, format::Format, response::Response};

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
    /// The response is an error.
    Response(response::Error),
    /// The data failed to decode.
    Decode(base64::DecodeError),
    /// The data URL is invalid.
    InvalidDataUrl,
}

impl error::Error for Error {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Url(e) => Some(e),
            Self::Request(e) => Some(e),
            Self::Response(e) => Some(e),
            Self::Decode(e) => Some(e),
            _ => None,
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
            Self::Decode(_) => f.write_str("decode error"),
            Self::InvalidDataUrl => f.write_str("invalid data URL"),
        }
    }
}
