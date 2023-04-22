#![warn(missing_docs)]

//! **noodles-fastq** handles the reading and writing of the FASTQ format.
//!
//! FASTQ is a text format with no formal specification and only has de facto rules. It typically
//! consists of a list of records, each with four lines: a definition (read name and description),
//! a sequence, a plus line, and quality scores.
//!
//! The read name is prefixed with an `@` (at sign) character and includes an optional description,
//! delimited by a space. The sequence is a list of bases encoded using IUPAC base symbols. The
//! plus line is effectively a separator, sometimes repeating the read name and optional
//! description, and is commonly discarded. The quality scores is list of Phred quality scores
//! offset by 33 and is parallel to each base in the sequence. That is, each record can be
//! described like the following:
//!
//! ```text
//! @<name>[ definition]
//! <sequence>
//! +[<name>[ definition]]
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
//!     .map(fastq::Reader::new)?;
//!
//! for result in reader.records() {
//!     let record = result?;
//!     // ...
//! }
//! # Ok::<(), io::Error>(())
//! ```

#[cfg(feature = "async")]
mod r#async;

pub mod fai;
mod indexer;
pub mod reader;
mod record;
mod writer;

pub use self::{indexer::Indexer, reader::Reader, record::Record, writer::Writer};

#[cfg(feature = "async")]
pub use self::r#async::{Reader as AsyncReader, Writer as AsyncWriter};

use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

/// Indexes a FASTQ file.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_fastq as fastq;
/// let index = fastq::index("sample.fastq")?;
/// # Ok::<(), io::Error>(())
/// ```
pub fn index<P>(src: P) -> io::Result<fai::Index>
where
    P: AsRef<Path>,
{
    let mut indexer = File::open(src).map(BufReader::new).map(Indexer::new)?;
    let mut index = Vec::new();

    while let Some(record) = indexer.index_record()? {
        index.push(record);
    }

    Ok(index)
}
