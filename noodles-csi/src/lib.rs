#![warn(missing_docs)]

//! **noodles-csi** handles the reading and writing of the coordinate-sorted index (CSI) format.

#[cfg(feature = "async")]
pub mod r#async;

pub mod binning_index;
pub mod io;
pub mod reader;
mod writer;

pub use self::{binning_index::BinningIndex, reader::Reader, writer::Writer};

use std::{fs::File, path::Path};

use self::binning_index::index::reference_sequence::index::BinnedIndex;

#[cfg(feature = "async")]
#[deprecated(
    since = "0.39.0",
    note = "Use `noodles_csi::r#async::io::Reader` instead."
)]
pub use self::r#async::io::Reader as AsyncReader;

#[cfg(feature = "async")]
#[deprecated(
    since = "0.39.0",
    note = "Use `noodles_csi::r#async::io::Writer` instead."
)]
pub use self::r#async::io::Writer as AsyncWriter;

/// A coordinate-sorted index (CSI).
pub type Index = binning_index::Index<BinnedIndex>;

/// Reads the entire contents of a coordinate-sorted index (CSI).
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
pub fn read<P>(src: P) -> std::io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    reader.read_index()
}

/// Writes a coordinate-sorted index (CSI) to a file.
///
/// This is a convenience function and is equivalent to creating a file at the given path and
/// writing the index.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_csi as csi;
/// let index = csi::Index::default();
/// csi::write("sample.bcf.csi", &index)?;
/// # Ok::<(), io::Error>(())
/// ```
pub fn write<P>(dst: P, index: &Index) -> std::io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::create(dst).map(Writer::new)?;
    writer.write_index(index)
}
