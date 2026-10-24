//! Binning index reference sequence index.

mod binned_index;
mod linear_index;

use noodles_bgzf as bgzf;
use noodles_core::Position;

pub use self::{binned_index::BinnedIndex, linear_index::LinearIndex};
use super::bin::Chunk;

/// A binning index reference sequence index.
pub trait Index {
    /// Returns the start virtual position of the first record in the bin that contains the given
    /// start position.
    fn min_offset(&self, min_shift: u8, depth: u8, start: Position) -> bgzf::VirtualPosition;

    /// Returns the start virtual position of the last first record.
    fn last_first_start_position(&self) -> Option<bgzf::VirtualPosition>;

    /// Finalizes record offsets once every record is added.
    ///
    /// `window_offsets` has, for each linear window (`1 << min_shift` positions), the start of the
    /// first record overlapping it, or `None`. The default implementation leaves the index
    /// unchanged.
    fn finish(&mut self, _depth: u8, _window_offsets: &[Option<bgzf::VirtualPosition>]) {}

    /// Adds a record to the index.
    fn update(&mut self, min_shift: u8, depth: u8, start: Position, end: Position, chunk: Chunk);
}
