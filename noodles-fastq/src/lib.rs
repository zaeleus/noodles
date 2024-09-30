#![warn(missing_docs)]

//! **noodles-fastq** handles the reading and writing of the FASTQ format.
//!
//! FASTQ is a text format with no formal specification and only has de facto rules. It typically
//! consists of a list of records, each with four lines: a definition (read name and description),
//! a sequence, a plus line, and quality scores.
//!
//! The read name is prefixed with an `@` (at sign) character and includes an optional description,
//! delimited by a space (` `) or horizontal tab (`\t`). The sequence is a list of bases encoded
//! using IUPAC base symbols. The plus line is effectively a separator, sometimes repeating the
//! read name and optional description, and is commonly discarded. The quality scores is list of
//! Phred quality scores (commonly but not guaranteed to be offset by 33) and is parallel to each
//! base in the sequence. That is, each record can be described like the following:
//!
//! ```text
//! @<name>[< |\t>description]
//! <sequence>
//! +[<name>[< |\t>description]]
//! <quality scores>
//! ```
//!
//! # Examples
//!
//! ## Read all records from a file
//!
//! ```no_run
//! # use std::{fs::File, io::{self, BufReader}};
//! use noodles_fastq as fastq;
//!
//! let mut reader = File::open("sample.fq")
//!     .map(BufReader::new)
//!     .map(fastq::io::Reader::new)?;
//!
//! for result in reader.records() {
//!     let record = result?;
//!     // ...
//! }
//! # Ok::<(), io::Error>(())
//! ```

#[cfg(feature = "async")]
pub mod r#async;

pub mod fai;
pub mod io;
pub mod record;

pub use self::record::Record;

#[deprecated(since = "0.11.0", note = "Use `noodles_fastq::io::index` instead.")]
pub use self::io::index;

#[deprecated(since = "0.11.0", note = "Use `noodles_fastq::io::Indexer` instead.")]
pub use self::io::Indexer;

#[deprecated(since = "0.11.0", note = "Use `noodles_fastq::io::reader` instead.")]
pub use self::io::reader;

#[deprecated(since = "0.11.0", note = "Use `noodles_fastq::io::Reader` instead.")]
pub use self::io::Reader;

//#[deprecated(since = "0.11.0", note = "Use `noodles_fastq::io::Writer` instead.")]
//pub use self::io::Writer;

#[cfg(feature = "async")]
#[deprecated(
    since = "0.15.0",
    note = "Use `noodles_fastq::r#async::io::Reader` instead."
)]
pub use self::r#async::io::Reader as AsyncReader;

#[deprecated(
    since = "0.15.0",
    note = "Use `noodles_fastq::r#async::io::Writer` instead."
)]
pub use self::r#async::io::Writer as AsyncWriter;
