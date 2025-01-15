//! CRAM index.

#[cfg(feature = "async")]
pub mod r#async;

pub mod io;
pub mod record;

pub use self::record::Record;

#[deprecated(since = "0.74.0", note = "Use `crai::io::Reader` instead.")]
pub use self::io::Reader;

#[deprecated(since = "0.74.0", note = "Use `crai::io::Writer` instead.")]
pub use self::io::Writer;

#[cfg(feature = "async")]
pub use self::r#async::{Reader as AsyncReader, Writer as AsyncWriter};

use std::{fs::File, path::Path};

/// A CRAM index.
pub type Index = Vec<Record>;

/// Reads the entire contents of a CRAM index.
///
/// This is a convenience function and is equivalent to opening the file at the given path and
/// reading the index.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_cram::crai;
/// let index = crai::read("sample.cram.crai")?;
/// # Ok::<(), io::Error>(())
/// ```
pub fn read<P>(src: P) -> std::io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    reader.read_index()
}

/// Writes a CRAM index to a file.
///
/// This is a convenience function and is equivalent to creating the file at the given path and
/// writing the index.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_cram::crai;
/// let index = crai::Index::default();
/// crai::write("sample.cram.crai", &index)?;
/// # Ok::<(), io::Error>(())
/// ```
pub fn write<P>(dst: P, index: &[Record]) -> std::io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::open(dst).map(Writer::new)?;
    writer.write_index(index)
}
