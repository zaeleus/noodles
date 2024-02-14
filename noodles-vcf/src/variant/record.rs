//! Variant record.

mod alternate_bases;
mod filters;
mod ids;
pub mod info;
pub mod samples;

pub use self::{
    alternate_bases::AlternateBases, filters::Filters, ids::Ids, info::Info, samples::Samples,
};
