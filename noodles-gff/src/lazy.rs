//! Lazily-evaluated GFF lines.

pub mod line;
pub mod record;

pub use self::{line::Line, record::Record};
