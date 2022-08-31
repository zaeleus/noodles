//! gzip index.
//!
//! A [gzip index] (GZI) is a list of compressed and uncompressed offset pairs for a gzipped file.
//!
//! [gzip index]: http://www.htslib.org/doc/bgzip.html#GZI_FORMAT

#[cfg(feature = "async")]
mod r#async;

mod reader;

pub use self::reader::Reader;

#[cfg(feature = "async")]
pub use self::r#async::Reader as AsyncReader;

use std::{fs::File, io, path::Path};

/// A gzip index.
pub type Index = Vec<(u64, u64)>;

/// Reads the entire contents of a GZ index.
///
/// This is a convenience function and is equivalent to opening the given path and reading the
/// index.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_bgzf::gzi;
/// let index = gzi::read("in.gz.gzi")?;
/// # Ok::<_, io::Error>(())
/// ```
pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    reader.read_index()
}
