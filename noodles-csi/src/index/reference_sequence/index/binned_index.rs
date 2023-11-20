use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use noodles_core::Position;

use super::Index;
use crate::index::reference_sequence::{bin::Chunk, parent_id, reg2bin};

/// A binned index.
pub type BinnedIndex = IndexMap<usize, bgzf::VirtualPosition>;

impl Index for BinnedIndex {
    fn min_offset(&self, min_shift: u8, depth: u8, start: Position) -> bgzf::VirtualPosition {
        let end = start;
        let mut bin_id = reg2bin(start, end, min_shift, depth);

        loop {
            if let Some(position) = self.get(&bin_id) {
                return *position;
            }

            bin_id = match parent_id(bin_id) {
                Some(id) => id,
                None => break,
            }
        }

        bgzf::VirtualPosition::default()
    }

    fn last_first_start_position(&self) -> Option<bgzf::VirtualPosition> {
        self.values().max().copied()
    }

    fn update(&mut self, min_shift: u8, depth: u8, start: Position, end: Position, chunk: Chunk) {
        let bin_id = reg2bin(start, end, min_shift, depth);

        self.entry(bin_id)
            .and_modify(|loffset| {
                if chunk.start() < *loffset {
                    *loffset = chunk.start();
                }
            })
            .or_insert(chunk.start());
    }
}
