//! BAM index.

mod format;
mod indexer;

use noodles_csi as csi;

pub use self::{format::Format, indexer::Indexer};
use crate::bai;

/// A BAM binning index.
pub enum Index {
    /// A BAM index.
    Bai(bai::Index),
    /// A coordinate-sorted index (CSI).
    Csi(csi::Index),
}
