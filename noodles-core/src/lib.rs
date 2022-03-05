#![warn(missing_docs)]

//! **noodles-core** contains shared structures and behavior among noodles libraries.

pub mod position;
pub mod region;

pub use self::{position::Position, region::Region};
