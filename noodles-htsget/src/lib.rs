#![warn(missing_docs)]

//! **noodles-htsget** is an htsget client.

pub(crate) mod chunks;
mod client;
mod format;
pub mod reads;
pub(crate) mod request;
mod response;
mod ticket;
pub mod variants;

pub(crate) use self::ticket::Ticket;
pub use self::{client::Client, format::Format, response::Response};

use std::{
    error, fmt,
    ops::{Bound, RangeBounds},
};

use noodles_core::Position;

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

impl error::Error for Error {}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Input => f.write_str("invalid input"),
            Self::Url(e) => write!(f, "URL error: {}", e),
            Self::Request(e) => write!(f, "request error: {}", e),
            Self::Response(e) => e.fmt(f),
            Self::Decode(e) => write!(f, "decode error: {}", e),
            Self::InvalidDataUrl => f.write_str("invalid data URL"),
        }
    }
}

// Resolves a 1-based [start, end] interval as a 0-based [start, end) interval.
fn resolve_interval<B>(interval: B) -> (Option<usize>, Option<usize>)
where
    B: RangeBounds<Position>,
{
    let start = match interval.start_bound() {
        Bound::Included(position) => Some(usize::from(*position) - 1),
        Bound::Excluded(position) => Some(usize::from(*position)),
        Bound::Unbounded => None,
    };

    let end = match interval.end_bound() {
        Bound::Included(position) => Some(usize::from(*position)),
        Bound::Excluded(position) => Some(usize::from(*position) - 1),
        Bound::Unbounded => None,
    };

    (start, end)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resolve_interval() -> std::result::Result<(), noodles_core::position::TryFromIntError> {
        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;

        // Range
        assert_eq!(resolve_interval(start..end), (Some(7), Some(12)));

        // RangeFrom
        assert_eq!(resolve_interval(start..), (Some(7), None));

        // RangeFull
        assert_eq!(resolve_interval(..), (None, None));

        // RangeInclusive
        assert_eq!(resolve_interval(start..=end), (Some(7), Some(13)));

        // RangeTo
        assert_eq!(resolve_interval(..end), (None, Some(12)));

        // RangeToInclusive
        assert_eq!(resolve_interval(..=end), (None, Some(13)));

        Ok(())
    }
}
