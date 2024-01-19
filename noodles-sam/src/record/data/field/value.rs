//! SAM record data field value and types.

pub mod base_modifications;
pub mod hex;

pub use self::{base_modifications::BaseModifications, hex::Hex};
