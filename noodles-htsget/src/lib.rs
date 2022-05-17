#![warn(missing_docs)]

//! **noodles-htsget** is an htsget client.

pub(crate) mod chunks;
mod client;
mod format;
pub mod reads;
mod response;
mod ticket;
pub mod variants;

pub(crate) use self::ticket::Ticket;
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
    /// The data failed to decode.
    Decode(base64::DecodeError),
    /// The data URL is invalid.
    InvalidDataUrl,
}

impl error::Error for Error {}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Input => f.write_str("invalid input"),
            Self::Url(e) => write!(f, "URL error: {}", e),
            Self::Request(e) => write!(f, "request error: {}", e),
            Self::Decode(e) => write!(f, "decode error: {}", e),
            Self::InvalidDataUrl => f.write_str("invalid data URL"),
        }
    }
}
