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

    /// Sets each bin's offset to that of the first record overlapping the bin's first window, as
    /// HTSlib does (`update_loff`).
    ///
    /// That record may be in an ancestor bin (e.g., a long record that started earlier), which the
    /// offsets from [`Self::update`] miss. A window that no record overlaps takes the offset of the
    /// previous one (or, before the first record, the first record's), which is still no later
    /// than any record at or after it.
    fn finish(&mut self, depth: u8, window_offsets: &[Option<bgzf::VirtualPosition>]) {
        let Some(first) = window_offsets.iter().find_map(|offset| *offset) else {
            return;
        };

        let mut prev = first;

        let filled: Vec<_> = window_offsets
            .iter()
            .map(|offset| {
                prev = offset.unwrap_or(prev);
                prev
            })
            .collect();

        for (&id, offset) in self.iter_mut() {
            // A bin always overlaps a record's window, so this only guards malformed input.
            *offset = filled
                .get(first_window(id, depth))
                .copied()
                .unwrap_or_default();
        }
    }

    fn update(&mut self, min_shift: u8, depth: u8, start: Position, end: Position, chunk: Chunk) {
        let mut bin_id = reg2bin(start, end, min_shift, depth);

        loop {
            self.entry(bin_id)
                .and_modify(|loffset| {
                    if chunk.start() < *loffset {
                        *loffset = chunk.start();
                    }
                })
                .or_insert(chunk.start());

            bin_id = match parent_id(bin_id) {
                Some(id) => id,
                None => break,
            };
        }
    }
}

/// Returns the index of the first linear window of a bin.
fn first_window(id: usize, depth: u8) -> usize {
    // Bins at level `l` start at ID `((1 << 3l) - 1) / 7`, and each spans `8^(depth - l)` windows.
    let mut level = 0;
    let mut level_start = 0;

    while level < depth && id >= level_start + (1 << (3 * u32::from(level))) {
        level_start += 1 << (3 * u32::from(level));
        level += 1;
    }

    (id - level_start) << (3 * u32::from(depth - level))
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
    fn test_first_window() {
        // depth 2: bin 0, then 8 bins of 8 windows, then 64 bins of 1 window. depth 5 (BAI-like):
        // level 4 starts at 585, and the leaves at 4681.
        for (id, depth, expected) in [
            (1, 2, 0),
            (2, 2, 8),
            (8, 2, 56),
            (9, 2, 0),
            (72, 2, 63),
            (586, 5, 8),
            (4681, 5, 0),
            (4781, 5, 100),
        ] {
            assert_eq!(first_window(id, depth), expected, "{id} {depth}");
        }
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
