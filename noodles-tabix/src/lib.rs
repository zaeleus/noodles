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
