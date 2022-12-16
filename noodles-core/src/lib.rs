#![warn(missing_docs)]

//! **noodles-core** contains shared structures and behavior among noodles libraries.

pub mod error;
pub mod position;
pub mod region;

pub use self::{error::Error, position::Position, region::Region};
