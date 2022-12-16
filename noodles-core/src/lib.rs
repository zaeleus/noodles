#![warn(missing_docs)]

//! **noodles-core** contains shared structures and behavior among noodles libraries.

pub mod error;
pub mod position;
pub mod region;

pub use self::{error::Error, position::Position, region::Region};

/// A specialized [std::result::Result] type for results in noodles.
pub type Result<T> = std::result::Result<T, error::Error>;
