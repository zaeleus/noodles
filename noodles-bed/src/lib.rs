#![warn(missing_docs)]

//! **noodles-bed** handles the reading and writing of the BED (Browser Extensible Data) format.

pub mod io;
pub mod record;
mod writer;

pub use self::{record::Record, writer::Writer};

#[deprecated(since = "0.14.0", note = "Use `noodles_bed::io::Reader` instead.")]
pub use self::io::Reader;
