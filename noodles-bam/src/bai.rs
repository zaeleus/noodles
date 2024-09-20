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
//! [querying]: crate::io::Reader::query
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

pub mod io;

#[deprecated(since = "0.68.0", note = "Use `bai::io::Reader` instead.")]
pub use self::io::Reader;

#[deprecated(since = "0.68.0", note = "Use `bai::io::Writer` instead.")]
pub use self::io::Writer;

#[cfg(feature = "async")]
#[deprecated(since = "0.68.0", note = "Use `bai::r#async::io::Reader` instead.")]
pub use self::r#async::Reader as AsyncReader;

#[cfg(feature = "async")]
pub use self::r#async::Writer as AsyncWriter;

use std::{
    fs::File,
    io::{BufReader, BufWriter},
    path::Path,
};

use noodles_csi::binning_index::{self, index::reference_sequence::index::LinearIndex};

const DEPTH: u8 = 5;

static MAGIC_NUMBER: &[u8] = b"BAI\x01";

/// A BAM index.
pub type Index = binning_index::Index<LinearIndex>;

/// Reads the entire contents of a BAM index.
///
/// This is a convenience function and is equivalent to opening the file at the given path, reading
/// the header, and reading the index.
///
/// # Examples
///
/// ```no_run
/// use noodles_bam::bai;
/// let index = bai::read("sample.bam.bai")?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn read<P>(src: P) -> std::io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(BufReader::new).map(Reader::new)?;
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
/// use noodles_bam::bai;
/// let index = bai::Index::default();
/// bai::write("sample.bam.bai", &index)?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn write<P>(dst: P, index: &Index) -> std::io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::create(dst).map(BufWriter::new).map(Writer::new)?;
    writer.write_index(index)
}
