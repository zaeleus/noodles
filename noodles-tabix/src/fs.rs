//! tabix filesystem operations.

use std::{fs::File, io, path::Path};

use super::{
    io::{Reader, Writer},
    Index,
};

/// Reads the entire contents of a tabix index.
///
/// This is a convenience function and is equivalent to opening the file at the given path and
/// reading the index.
///
/// # Examples
///
/// ```no_run
/// use noodles_tabix as tabix;
/// let index = tabix::fs::read("sample.vcf.gz.tbi")?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    reader.read_index()
}

/// Writes a tabix index to a file.
///
/// This is a convenience function and is equivalent to creating a file at the given path and
/// writing the index.
///
/// # Examples
///
/// ```no_run
/// use noodles_csi::binning_index::index::Header;
/// use noodles_tabix as tabix;
/// let index = tabix::Index::builder().set_header(Header::default()).build();
/// tabix::fs::write("sample.vcf.gz.tbi", &index)?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn write<P>(dst: P, index: &Index) -> io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::create(dst).map(Writer::new)?;
    writer.write_index(index)?;
    Ok(())
}
