#![warn(missing_docs)]

//! **noodles-gtf** handles the reading and writing of the Gene Transfer Format (GTF).

pub mod io;
pub mod line;
pub mod record;
mod writer;

pub use self::{line::Line, record::Record, writer::Writer};

#[deprecated(since = "0.32.0", note = "Use `noodles_gtf::io::Reader` instead.")]
pub use self::io::Reader;
