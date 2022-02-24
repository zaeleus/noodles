//! FASTA index (FAI) and fields.

#[cfg(feature = "async")]
mod r#async;

mod reader;
mod record;
mod writer;

pub use self::{reader::Reader, record::Record, writer::Writer};

#[cfg(feature = "async")]
pub use self::r#async::Reader as AsyncReader;

use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

/// A FASTA index.
pub type Index = Vec<Record>;

/// Reads the entire contents of a FASTA index.
///
/// This is a convenience function and is equivalent to opening the file at the given path and
/// parsing each record.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_fasta::fai;
/// let index = fai::read("reference.fa.fai")?;
/// # Ok::<(), io::Error>(())
/// ```
pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(BufReader::new).map(Reader::new)?;
    reader.read_index()
}
