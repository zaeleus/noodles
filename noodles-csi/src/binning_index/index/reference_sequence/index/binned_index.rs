use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use noodles_core::Position;

use super::Index;
use crate::binning_index::index::reference_sequence::{bin::Chunk, parent_id, reg2bin};

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_min_offset() -> Result<(), noodles_core::position::TryFromIntError> {
        const MIN_SHIFT: u8 = 4;
        const DEPTH: u8 = 2;

        let index: BinnedIndex = [(1, bgzf::VirtualPosition::from(233))]
            .into_iter()
            .collect();

        let start = Position::try_from(8)?;
        assert_eq!(
            index.min_offset(MIN_SHIFT, DEPTH, start),
            bgzf::VirtualPosition::from(233)
        );

        let start = Position::try_from(144)?;
        assert_eq!(
            index.min_offset(MIN_SHIFT, DEPTH, start),
            bgzf::VirtualPosition::default()
        );

        Ok(())
    }

    #[test]
    fn test_last_first_start_position() {
        let index: BinnedIndex = [
            (0, bgzf::VirtualPosition::from(8)),
            (2, bgzf::VirtualPosition::from(21)),
            (9, bgzf::VirtualPosition::from(13)),
        ]
        .into_iter()
        .collect();

        assert_eq!(
            index.last_first_start_position(),
            Some(bgzf::VirtualPosition::from(21))
        );
    }
}
