//! BAM index (BAI) and fields.
//!
//! A BAM index (BAI) is used with an associated coordinate-sorted BAM file that allows random
//! access to records, e.g., [querying].
//!
//! The index contains a list of reference sequences parallel to the one defined in the BAM file.
//! Each indexed reference sequence has a calculated set of hierarchical bins at different
//! granularities. The bins then define a list of physical file positions in the BAM to search for
//! overlapping records.
//!
//! When reading entire BAM files sequentially, a BAM index is not necessary.
//!
//! [querying]: crate::Reader::query
//!
//! # Examples
//!
//! ## Reading a BAM index
//!
//! ```no_run
//! # use std::io;
//! use noodles_bam::bai;
//! let index = bai::read("sample.bam.bai")?;
//! # Ok::<(), io::Error>(())
//! ```

#[cfg(feature = "async")]
pub mod r#async;

pub mod index;
mod reader;
mod writer;

pub use self::{index::Index, reader::Reader, writer::Writer};

#[cfg(feature = "async")]
pub use self::r#async::{Reader as AsyncReader, Writer as AsyncWriter};

use std::{fs::File, io, path::Path};

use self::index::reference_sequence::Bin;

static MAGIC_NUMBER: &[u8] = b"BAI\x01";

/// Reads the entire contents of a BAM index.
///
/// This is a convenience function and is equivalent to opening the file at the given path, reading
/// the header, and reading the index.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_bam::bai;
/// let index = bai::read("sample.bam.bai")?;
/// # Ok::<(), io::Error>(())
/// ```
pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    reader.read_header()?;
    reader.read_index()
}

/// Writes a BAM index to a file.
///
/// This is a convenience function and is equivalent to creating a file at the given path, writing
/// the header, and writing the index.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_bam::bai;
/// let index = bai::Index::default();
/// bai::write("sample.bam.bai", &index)?;
/// # Ok::<(), io::Error>(())
/// ```
pub fn write<P>(dst: P, index: &Index) -> io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::create(dst).map(Writer::new)?;
    writer.write_header()?;
    writer.write_index(index)
}
