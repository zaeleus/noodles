#![warn(missing_docs)]

//! **noodles-csi** handles the reading and writing of the coordinate-sorted index (CSI) format.

#[cfg(feature = "async")]
pub mod r#async;

pub mod binning_index;
pub mod index;
mod reader;
mod writer;

pub use self::{binning_index::BinningIndex, index::Index, reader::Reader, writer::Writer};

#[deprecated(
    since = "0.4.0",
    note = "Use `csi::binning_index::ReferenceSequenceExt` instead."
)]
pub use binning_index::ReferenceSequenceExt as BinningIndexReferenceSequence;

#[cfg(feature = "async")]
pub use self::r#async::{Reader as AsyncReader, Writer as AsyncWriter};

use std::{fs::File, io, path::Path};

static MAGIC_NUMBER: &[u8] = b"CSI\x01";

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
pub fn read<P>(src: P) -> io::Result<Index>
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
pub fn write<P>(dst: P, index: &Index) -> io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::create(dst).map(Writer::new)?;
    writer.write_index(index)
}
