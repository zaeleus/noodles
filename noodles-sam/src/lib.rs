//! **noodles-sam** handles the reading and writing of the SAM (Sequence Alignment/Map) format.
//!
//! SAM is a format typically used to store biological sequences, either mapped to a reference
//! sequence or unmapped. It has two sections: a header and a list of records.
//!
//! The header mostly holds meta information about the data: a header describing the file
//! format version, reference sequences reads map to, read groups reads belong to, programs that
//! previously manipulated the data, and free-form comments. The header is optional and may be
//! empty.
//!
//! Each record represents a read, a linear alignment of a segment. Records have fields describing
//! how a read was mapped (or not) to a reference sequence.
//!
//! # Examples
//!
//! ## Read all records from a file
//!
//! ```no_run
//! use noodles_sam as sam;
//!
//! let mut reader = sam::io::reader::Builder::default().build_from_path("sample.sam")?;
//! let header = reader.read_header()?;
//!
//! for result in reader.records() {
//!     let record = result?;
//!     // ...
//! }
//! # Ok::<_, std::io::Error>(())
//! ```

#[cfg(feature = "async")]
pub mod r#async;

pub mod alignment;
pub mod fs;
pub mod header;
pub mod io;
pub mod record;

pub use self::{header::Header, record::Record};
