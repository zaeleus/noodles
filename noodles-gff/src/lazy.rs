//! Lazily-evaluated GFF lines.

mod directive;
pub mod line;
pub mod record;

pub use self::{directive::Directive, line::Line, record::Record};
