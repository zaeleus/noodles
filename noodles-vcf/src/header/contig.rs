//! VCF header contig record components.

pub mod name;
mod tag;

pub use self::{name::Name, tag::Tag};
