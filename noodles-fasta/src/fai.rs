//! FASTA index (FAI) and fields.

#[cfg(feature = "async")]
mod r#async;

mod index;
pub mod io;
mod record;

pub use self::{index::Index, record::Record};

#[deprecated(since = "0.44.0", note = "Use `fai::io::Reader` instead.")]
pub use self::io::Reader;

#[deprecated(since = "0.44.0", note = "Use `fai::io::Writer` instead.")]
pub use self::io::Writer;

#[cfg(feature = "async")]
pub use self::r#async::Reader as AsyncReader;

use std::{fs::File, io::BufReader, path::Path};

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
pub fn read<P>(src: P) -> std::io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(BufReader::new).map(Reader::new)?;
    reader.read_index()
}
