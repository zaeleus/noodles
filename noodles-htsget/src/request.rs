//! htsget request.

mod builder;
mod class;
mod kind;
mod payload;

pub use self::class::Class;

pub(crate) use self::{builder::Builder, kind::Kind, payload::Payload};
