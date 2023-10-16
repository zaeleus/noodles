//! Lazily-evaluated GFF lines.

mod line;
pub(crate) mod record;

pub use self::{line::Line, record::Record};
