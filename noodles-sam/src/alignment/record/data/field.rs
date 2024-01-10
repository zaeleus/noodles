//! Alignment record data field.

pub mod tag;
mod ty;
pub mod value;

pub use self::{tag::Tag, ty::Type, value::Value};
