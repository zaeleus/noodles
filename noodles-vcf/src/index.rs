//! VCF index.

mod format;
pub mod indexer;

use std::io;

use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::{
    self as csi, BinningIndex, binning_index::index::reference_sequence::bin::Chunk,
};
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

impl BinningIndex for Index {
    fn min_shift(&self) -> u8 {
        match self {
            Self::Csi(index) => index.min_shift(),
            Self::Tabix(index) => index.min_shift(),
        }
    }

    fn depth(&self) -> u8 {
        match self {
            Self::Csi(index) => index.depth(),
            Self::Tabix(index) => index.depth(),
        }
    }

    fn header(&self) -> Option<&csi::binning_index::index::Header> {
        match self {
            Self::Csi(index) => index.header(),
            Self::Tabix(index) => index.header(),
        }
    }

    fn reference_sequences(
        &self,
    ) -> Box<dyn Iterator<Item = &dyn csi::binning_index::ReferenceSequence> + '_> {
        match self {
            Self::Csi(index) => BinningIndex::reference_sequences(index),
            Self::Tabix(index) => BinningIndex::reference_sequences(index),
        }
    }

    fn unplaced_unmapped_record_count(&self) -> Option<u64> {
        match self {
            Self::Csi(index) => index.unplaced_unmapped_record_count(),
            Self::Tabix(index) => index.unplaced_unmapped_record_count(),
        }
    }

    fn query(&self, reference_sequence_id: usize, interval: Interval) -> io::Result<Vec<Chunk>> {
        match self {
            Self::Csi(index) => index.query(reference_sequence_id, interval),
            Self::Tabix(index) => index.query(reference_sequence_id, interval),
        }
    }

    fn last_first_record_start_position(&self) -> Option<bgzf::VirtualPosition> {
        match self {
            Self::Csi(index) => index.last_first_record_start_position(),
            Self::Tabix(index) => index.last_first_record_start_position(),
        }
    }
}
