#![warn(missing_docs)]

//! **noodles-tabix** handles the reading and writing of the [tabix format].
//!
//! A tabix (TBI) is an index file typically used to allow random access of an accompanied file
//! that is
//!
//!   1) bgzipped,
//!   2) tab-delimited,
//!   3) grouped by reference sequence name, and
//!   4) coordinate sorted by start position.
//!
//! It can be used to find relevant records for a given genomic region.
//!
//! [tabix format]: https://samtools.github.io/hts-specs/tabix.pdf
//!
//! # Examples
//!
//! ## Read a tabix file
//!
//! ```no_run
//! # use std::io;
//! use noodles_tabix as tabix;
//! let index = tabix::read("sample.vcf.gz.tbi")?;
//! # Ok::<(), io::Error>(())
//! ```

pub mod index;
mod reader;
mod writer;

pub use self::{index::Index, reader::Reader, writer::Writer};

use std::{fs::File, io, path::Path};

static MAGIC_NUMBER: &[u8] = b"TBI\x01";

/// Reads the entire contents of a tabix index.
///
/// This is a convenience function and is equivalent to opening the file at the given path and
/// reading the index.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_tabix as tabix;
/// let index = tabix::read("sample.vcf.gz.tbi")?;
/// # Ok::<(), io::Error>(())
/// ```
pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    reader.read_index()
}
