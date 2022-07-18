//! VCF header genotype format record components.

pub mod key;
mod tag;
pub mod ty;

pub use self::{key::Key, tag::Tag, ty::Type};
