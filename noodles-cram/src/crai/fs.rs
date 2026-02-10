//! CRAM index filesystem operations.

use std::{fs::File, path::Path};

use super::{
    Index, Record,
    io::{Reader, Writer},
};

/// Reads the entire contents of a CRAM index.
///
/// This is a convenience function and is equivalent to opening the file at the given path and
/// reading the index.
///
/// # Examples
///
/// ```no_run
/// use noodles_cram::crai;
/// let index = crai::fs::read("sample.cram.crai")?;
/// # Ok::<(), std::io::Error>(())
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
/// use noodles_cram::crai;
/// let index = crai::Index::default();
/// crai::fs::write("sample.cram.crai", &index)?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn write<P>(dst: P, index: &[Record]) -> std::io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::create(dst).map(Writer::new)?;
    writer.write_index(index)
}
