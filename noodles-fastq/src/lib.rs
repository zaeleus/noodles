#![deny(missing_docs)]

//! **noodles-fastq** handles the reading and writing of the FASTQ format.
//!
//! FASTQ is a text format with no formal specification and only has de facto rules. It typically
//! consists of a list of records, each with four lines: a read name, a sequence, a plus line, and
//! quality scores.
//!
//! The read name is prefixed with an `@` (at sign) character. The sequence is a list of bases
//! encoded using IUPAC base symbols. The plus line is effectively a separator, sometimes repeating
//! the read name, and is commonly discarded. The quality scores is list of Phred quality scores
//! offset by 33 and is parallel to a base in the sequence.
//!
//! # Examples
//!
//! ## Read all records from a file
//!
//! ```no_run
//! # use std::{fs::File, io::{self, BufReader}};
//! use noodles_fastq as fastq;
//!
//! let mut reader = File::open("sample.fastq").map(BufReader::new).map(fastq::Reader::new)?;
//!
//! for result in reader.records() {
//!     let record = result?;
//!     println!("{:?}", record);
//! }
//! # Ok::<(), io::Error>(())
//! ```

pub mod fai;
mod reader;
mod record;
mod writer;

pub use self::{reader::Reader, record::Record, writer::Writer};
