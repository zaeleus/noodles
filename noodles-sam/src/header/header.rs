//! SAM header header and fields.
//!
//! The namespace of this module is intentionally awkward to disambiguate a SAM header
//! ([`crate::Header`]) and a header record ([`crate::header::header::Header`]).

pub mod group_order;
pub mod sort_order;
pub mod subsort_order;
pub mod version;

pub use self::{
    group_order::GroupOrder, sort_order::SortOrder, subsort_order::SubsortOrder, version::Version,
};
