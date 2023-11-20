//! CSI reference sequence index.

mod binned_index;
mod linear_index;

use noodles_bgzf as bgzf;
use noodles_core::Position;

pub use self::{binned_index::BinnedIndex, linear_index::LinearIndex};
use super::bin::Chunk;

/// A CSI reference sequence index.
pub trait Index {
    /// Returns the start virtual position of the first record in the bin that contains the given
    /// start position.
    fn min_offset(&self, min_shift: u8, depth: u8, start: Position) -> bgzf::VirtualPosition;

    /// Returns the start virtual position of the last first record.
    fn last_first_start_position(&self) -> Option<bgzf::VirtualPosition>;

    /// Adds a record to the index.
    fn update(&mut self, min_shift: u8, depth: u8, start: Position, end: Position, chunk: Chunk);
}
