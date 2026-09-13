//! BAM index.

mod format;
pub mod indexer;

use std::io;

use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::{
    self as csi, BinningIndex, binning_index::index::reference_sequence::bin::Chunk,
};

pub use self::{format::Format, indexer::Indexer};
use crate::bai;

/// A BAM binning index.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Index {
    /// A BAM index.
    Bai(bai::Index),
    /// A coordinate-sorted index (CSI).
    Csi(csi::Index),
}

impl BinningIndex for Index {
    fn min_shift(&self) -> u8 {
        match self {
            Self::Bai(index) => index.min_shift(),
            Self::Csi(index) => index.min_shift(),
        }
    }

    fn depth(&self) -> u8 {
        match self {
            Self::Bai(index) => index.depth(),
            Self::Csi(index) => index.depth(),
        }
    }

    fn header(&self) -> Option<&csi::binning_index::index::Header> {
        match self {
            Self::Bai(index) => index.header(),
            Self::Csi(index) => index.header(),
        }
    }

    fn reference_sequences(
        &self,
    ) -> Box<dyn Iterator<Item = &dyn csi::binning_index::ReferenceSequence> + '_> {
        match self {
            Self::Bai(index) => BinningIndex::reference_sequences(index),
            Self::Csi(index) => BinningIndex::reference_sequences(index),
        }
    }

    fn unplaced_unmapped_record_count(&self) -> Option<u64> {
        match self {
            Self::Bai(index) => index.unplaced_unmapped_record_count(),
            Self::Csi(index) => index.unplaced_unmapped_record_count(),
        }
    }

    fn query(&self, reference_sequence_id: usize, interval: Interval) -> io::Result<Vec<Chunk>> {
        match self {
            Self::Bai(index) => index.query(reference_sequence_id, interval),
            Self::Csi(index) => index.query(reference_sequence_id, interval),
        }
    }

    fn last_first_record_start_position(&self) -> Option<bgzf::VirtualPosition> {
        match self {
            Self::Bai(index) => index.last_first_record_start_position(),
            Self::Csi(index) => index.last_first_record_start_position(),
        }
    }
}
