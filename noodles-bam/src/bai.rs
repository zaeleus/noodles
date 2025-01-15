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
//! let index = bai::fs::read("sample.bam.bai")?;
//! # Ok::<(), io::Error>(())
//! ```

#[cfg(feature = "async")]
pub mod r#async;

pub mod fs;
pub mod io;

#[deprecated(since = "0.73.0", note = "Use `bai::fs::read` instead.")]
pub use self::fs::read;

#[deprecated(since = "0.73.0", note = "Use `bai::fs::write` instead.")]
pub use self::fs::write;

#[deprecated(since = "0.68.0", note = "Use `bai::io::Reader` instead.")]
pub use self::io::Reader;

#[deprecated(since = "0.68.0", note = "Use `bai::io::Writer` instead.")]
pub use self::io::Writer;

#[cfg(feature = "async")]
#[deprecated(since = "0.68.0", note = "Use `bai::r#async::io::Reader` instead.")]
pub use self::r#async::Reader as AsyncReader;

#[cfg(feature = "async")]
#[deprecated(since = "0.68.0", note = "Use `bai::r#async::io::Writer` instead.")]
pub use self::r#async::Writer as AsyncWriter;

use noodles_csi::binning_index::{self, index::reference_sequence::index::LinearIndex};

const DEPTH: u8 = 5;

const MAGIC_NUMBER: [u8; 4] = *b"BAI\x01";

/// A BAM index.
pub type Index = binning_index::Index<LinearIndex>;
