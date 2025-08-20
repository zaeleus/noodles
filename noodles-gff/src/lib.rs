//! **noodles-gff** handles the reading and writing of the [GFF3 format][gff3-spec].
//!
//! GFF (Generic Feature Format) is a text-based format used to represent genomic features.
//!
//! [gff3-spec]: https://github.com/The-Sequence-Ontology/Specifications/blob/be6e1af7243ba4235c30b69660e2669e444e2f3e/gff3.md
//!
//! # Examples
//!
//! ## Read all records
//!
//! ```no_run
//! # use std::{fs::File, io::{self, BufReader}};
//! use noodles_gff as gff;
//!
//! let mut reader = File::open("annotations.gff3")
//!     .map(BufReader::new)
//!     .map(gff::io::Reader::new)?;
//!
//! for result in reader.record_bufs() {
//!     let record = result?;
//!     // ...
//! }
//! # Ok::<(), io::Error>(())
//! ```

#[cfg(feature = "async")]
pub mod r#async;

mod directive;
pub mod directive_buf;
pub mod feature;
pub mod io;
pub mod line;
pub mod line_buf;
pub mod record;

pub use self::{
    directive::Directive, directive_buf::DirectiveBuf, line::Line, line_buf::LineBuf,
    record::Record,
};
