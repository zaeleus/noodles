//! VCF index.

mod format;
mod indexer;

use noodles_csi as csi;
use noodles_tabix as tabix;

pub use self::{format::Format, indexer::Indexer};

/// A VCF binning index.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Index {
    /// A coordinate-sorted index (CSI).
    Csi(csi::Index),
    /// A tabix (TBI).
    Tabix(tabix::Index),
}
