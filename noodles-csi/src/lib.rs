#![warn(missing_docs)]

//! **noodles-csi** handles the reading and writing of the coordinate-sorted index (CSI) format.

pub mod index;
mod reader;
mod writer;

pub use self::{index::Index, reader::Reader, writer::Writer};

use std::{fs::File, io, path::Path};

static MAGIC_NUMBER: &[u8] = b"CSI\x01";

/// Reads the entire contents of an coordinate-sorted index (CSI).
///
/// This is a convenience function and is equivalent to opening the file at the given path and
/// reading the index.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_csi as csi;
/// let index = csi::read("sample.bcf.csi")?;
/// # Ok::<(), io::Error>(())
/// ```
pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    reader.read_index()
}
