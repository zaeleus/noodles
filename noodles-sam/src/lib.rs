#![warn(missing_docs)]

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
//! let mut reader = sam::reader::Builder::default().build_from_path("sample.sam")?;
//! let header = reader.read_header()?;
//!
//! for result in reader.records(&header) {
//!     let record = result?;
//!     // ...
//! }
//! # Ok::<_, std::io::Error>(())
//! ```

#[cfg(feature = "async")]
mod r#async;

pub mod alignment;
mod alignment_writer;
pub mod header;
pub mod indexed_reader;
pub mod io;
pub mod lazy;
pub mod reader;
pub mod record;
pub mod writer;

pub use self::{
    alignment_writer::AlignmentWriter, header::Header, indexed_reader::IndexedReader,
    reader::Reader, writer::Writer,
};

#[cfg(feature = "async")]
pub use self::r#async::{Reader as AsyncReader, Writer as AsyncWriter};
