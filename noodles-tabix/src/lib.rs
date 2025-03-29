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
//! use noodles_tabix as tabix;
//! let index = tabix::fs::read("sample.vcf.gz.tbi")?;
//! # Ok::<(), std::io::Error>(())
//! ```

#[cfg(feature = "async")]
pub mod r#async;

pub mod fs;
pub mod index;
pub mod io;

#[deprecated(since = "0.48.0", note = "Use `tabix::fs::read` instead.")]
pub use self::fs::read;

#[deprecated(since = "0.48.0", note = "Use `tabix::fs::write` instead.")]
pub use self::fs::write;

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
#[deprecated(
    since = "0.45.0",
    note = "Use `noodles_tabix::r#async::io::Writer` instead."
)]
pub use self::r#async::Writer as AsyncWriter;

use noodles_csi::binning_index::{self, index::reference_sequence::index::LinearIndex};

const MAGIC_NUMBER: [u8; 4] = *b"TBI\x01";

/// A tabix index.
pub type Index = binning_index::Index<LinearIndex>;
