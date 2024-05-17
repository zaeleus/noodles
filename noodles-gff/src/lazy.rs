//! Lazily-evaluated GFF lines.

mod line;
pub mod record;

pub use self::{line::Line, record::Record};
