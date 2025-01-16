//! gzip index filesystem operations.

use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

use super::{io::Reader, Index};

/// Reads the entire contents of a GZ index.
///
/// This is a convenience function and is equivalent to opening the given path and reading the
/// index.
///
/// # Examples
///
/// ```no_run
/// use noodles_bgzf::gzi;
/// let index = gzi::fs::read("in.gz.gzi")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(BufReader::new).map(Reader::new)?;
    reader.read_index()
}
