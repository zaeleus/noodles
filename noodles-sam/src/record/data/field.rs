//! SAM record data field and components.

pub mod tag;
pub mod ty;
pub mod value;

pub use self::{tag::Tag, ty::Type, value::Value};
