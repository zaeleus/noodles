#![warn(missing_docs)]

//! **noodles-csi** handles the reading and writing of the coordinate-sorted index (CSI) format.

#[cfg(feature = "async")]
pub mod r#async;

pub mod binning_index;
pub mod fs;
pub mod io;

use self::binning_index::index::reference_sequence::index::BinnedIndex;
pub use self::binning_index::BinningIndex;

#[deprecated(since = "0.42.0", note = "Use `csi::fs::read` instead.")]
pub use self::fs::read;

#[deprecated(since = "0.42.0", note = "Use `csi::fs::write` instead.")]
pub use self::fs::write;

#[deprecated(since = "0.39.0", note = "Use `noodles_csi::io::Reader` instead.")]
pub use self::io::Reader;

#[deprecated(since = "0.39.0", note = "Use `noodles_csi::io::Writer` instead.")]
pub use self::io::Writer;

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
