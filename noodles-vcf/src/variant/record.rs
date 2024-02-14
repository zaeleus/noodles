//! Variant record.

mod alternate_bases;
mod ids;
pub mod info;
pub mod samples;

pub use self::{alternate_bases::AlternateBases, ids::Ids, info::Info, samples::Samples};
