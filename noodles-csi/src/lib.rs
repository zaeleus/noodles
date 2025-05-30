//! **noodles-csi** handles the reading and writing of the coordinate-sorted index (CSI) format.

#[cfg(feature = "async")]
pub mod r#async;

pub mod binning_index;
pub mod fs;
pub mod io;

pub use self::binning_index::BinningIndex;
use self::binning_index::index::reference_sequence::index::BinnedIndex;

/// A coordinate-sorted index (CSI).
pub type Index = binning_index::Index<BinnedIndex>;
