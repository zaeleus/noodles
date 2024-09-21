#![warn(missing_docs)]

//! **noodles-tabix** handles the reading and writing of the [tabix format].
//!
//! A tabix (TBI) is an index file typically used to allow random access of an accompanied file
//! that is
//!
//!   1) bgzipped,
//!   2) tab-delimited,
//!   3) grouped by reference sequence name, and
//!   4) coordinate sorted by start position.
//!
//! It can be used to find relevant records for a given genomic region.
//!
//! [tabix format]: https://samtools.github.io/hts-specs/tabix.pdf
//!
//! # Examples
//!
//! ## Read a tabix file
//!
//! ```no_run
//! # use std::io;
//! use noodles_tabix as tabix;
//! let index = tabix::read("sample.vcf.gz.tbi")?;
//! # Ok::<(), io::Error>(())
//! ```

#[cfg(feature = "async")]
pub mod r#async;

pub mod index;
pub mod io;

#[deprecated(since = "0.45.0", note = "Use `noodles_tabix::io::Reader` instead.")]
pub use self::io::Reader;

#[deprecated(since = "0.45.0", note = "Use `noodles_tabix::io::Writer` instead.")]
pub use self::io::Writer;

#[cfg(feature = "async")]
#[deprecated(
    since = "0.45.0",
    note = "Use `noodles_tabix::r#async::io::Reader` instead."
)]
pub use self::r#async::Reader as AsyncReader;

#[cfg(feature = "async")]
pub use self::r#async::Writer as AsyncWriter;

use std::{fs::File, path::Path};

use noodles_csi::binning_index::{self, index::reference_sequence::index::LinearIndex};

static MAGIC_NUMBER: &[u8] = b"TBI\x01";

/// A tabix index.
pub type Index = binning_index::Index<LinearIndex>;

/// Reads the entire contents of a tabix index.
///
/// This is a convenience function and is equivalent to opening the file at the given path and
/// reading the index.
///
/// # Examples
///
/// ```no_run
/// use noodles_tabix as tabix;
/// let index = tabix::read("sample.vcf.gz.tbi")?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn read<P>(src: P) -> std::io::Result<Index>
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
///
/// let index = tabix::Index::builder().set_header(Header::default()).build();
/// tabix::write("sample.vcf.gz.tbi", &index)?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn write<P>(dst: P, index: &Index) -> std::io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::create(dst).map(Writer::new)?;
    writer.write_index(index)?;
    Ok(())
}
