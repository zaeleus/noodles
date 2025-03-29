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
//! let mut reader = bam::io::reader::Builder::default().build_from_path("sample.bam")?;
//! let header = reader.read_header()?;
//!
//! for result in reader.records() {
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
//! let mut reader = bam::io::indexed_reader::Builder::default().build_from_path("sample.bam")?;
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
pub mod r#async;

pub mod bai;
pub mod fs;
pub mod io;
pub mod record;

pub use self::record::Record;

#[cfg(feature = "async")]
#[deprecated(since = "0.68.0", note = "Use `bam::r#async::io::Reader` instead.")]
pub use self::r#async::io::Reader as AsyncReader;

#[cfg(feature = "async")]
#[deprecated(since = "0.68.0", note = "Use `bam::r#async::io::Writer` instead.")]
pub use self::r#async::io::Writer as AsyncWriter;

const MAGIC_NUMBER: [u8; 4] = *b"BAM\x01";
