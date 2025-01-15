//! FAI filesystem operations.

use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

use super::{io::Reader, Index};

/// Reads the entire contents of a FASTA index.
///
/// This is a convenience function and is equivalent to opening the file at the given path and
/// parsing each record.
///
/// # Examples
///
/// ```no_run
/// use noodles_fasta::fai;
/// let index = fai::fs::read("reference.fa.fai")?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(BufReader::new).map(Reader::new)?;
    reader.read_index()
}
