#![warn(missing_docs)]

//! **noodles-bam** handles the reading and writing of the BAM (Binary Alignment/Map) file format.
//!
//! The BAM format contains the same information as SAM (Sequence Alignment/Map), namely a SAM
//! header and a list of records.
//!
//! # Examples
//!
//! ## Read all records
//!
//! ```no_run
//! # use std::io;
//! use noodles_bam as bam;
//!
//! let mut reader = bam::reader::Builder::default().build_from_path("sample.bam")?;
//! let header = reader.read_header()?;
//!
//! for result in reader.records(&header) {
//!     let record = result?;
//!     // ...
//! }
//! # Ok::<_, io::Error>(())
//! ```
//!
//! ## Query records
//!
//! Querying allows filtering records by region. It requires an associated BAM index (BAI).
//!
//! ```no_run
//! # use std::fs::File;
//! use noodles_bam as bam;
//!
//! let mut reader = bam::indexed_reader::Builder::default().build_from_path("sample.bam")?;
//! let header = reader.read_header()?;
//!
//! let region = "sq0:5-8".parse()?;
//! let query = reader.query(&header, &region)?;
//!
//! for result in query {
//!     let record = result?;
//!     // ...
//! }
//! # Ok::<_, Box<dyn std::error::Error>>(())
//! ```

#[cfg(feature = "async")]
mod r#async;

pub mod bai;
pub mod indexed_reader;
pub mod lazy;
pub mod reader;
pub mod writer;

#[doc(hidden)]
pub mod record;

pub use self::{indexed_reader::IndexedReader, reader::Reader, writer::Writer};

#[cfg(feature = "async")]
pub use self::r#async::{Reader as AsyncReader, Writer as AsyncWriter};

static MAGIC_NUMBER: &[u8] = b"BAM\x01";
