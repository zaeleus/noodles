//! Binning index utilities.

mod reference_sequence_ext;

pub use reference_sequence_ext::ReferenceSequenceExt;

use std::io;

use noodles_bgzf as bgzf;
use noodles_core::region::Interval;

use super::index::reference_sequence::bin::Chunk;

/// A binning index.
pub trait BinningIndex {
    /// The returned output indexed reference sequence.
    type ReferenceSequence: ReferenceSequenceExt;

    /// Returns a list of indexed reference sequences.
    fn reference_sequences(&self) -> &[Self::ReferenceSequence];

    /// Returns the number of unplaced, unmapped records in the associated file.
    fn unplaced_unmapped_record_count(&self) -> Option<u64>;

    /// Returns the chunks that overlap with the given region.
    fn query<I>(&self, reference_sequence_id: usize, interval: I) -> io::Result<Vec<Chunk>>
    where
        I: Into<Interval>;

    /// Returns the start position of the first record in the last linear bin.
    ///
    /// This is the closest position to the unplaced, unmapped records, if any, that is available
    /// in an index.
    fn first_record_in_last_linear_bin_start_position(&self) -> Option<bgzf::VirtualPosition> {
        self.reference_sequences()
            .iter()
            .rev()
            .find_map(|rs| rs.first_record_in_last_linear_bin_start_position())
    }
}

/// Merges a list of chunks into a list of non-overlapping chunks.
///
/// This is the same as calling [`optimize_chunks`] with a `min_offset` of 0.
///
/// # Examples
///
/// ```
/// use noodles_bgzf as bgzf;
/// use noodles_csi::{
///     binning_index::merge_chunks,
///     index::reference_sequence::bin::Chunk,
/// };
///
/// let chunks = [
///     Chunk::new(bgzf::VirtualPosition::from(2), bgzf::VirtualPosition::from(3)),
///     Chunk::new(bgzf::VirtualPosition::from(5), bgzf::VirtualPosition::from(8)),
///     Chunk::new(bgzf::VirtualPosition::from(7), bgzf::VirtualPosition::from(13)),
///     Chunk::new(bgzf::VirtualPosition::from(21), bgzf::VirtualPosition::from(34)),
/// ];
///
/// let actual = merge_chunks(&chunks);
///
/// let expected = [
///     Chunk::new(bgzf::VirtualPosition::from(2), bgzf::VirtualPosition::from(3)),
///     Chunk::new(bgzf::VirtualPosition::from(5), bgzf::VirtualPosition::from(13)),
///     Chunk::new(bgzf::VirtualPosition::from(21), bgzf::VirtualPosition::from(34)),
/// ];
///
/// assert_eq!(actual, expected);
/// ```
pub fn merge_chunks(chunks: &[Chunk]) -> Vec<Chunk> {
    optimize_chunks(chunks, bgzf::VirtualPosition::default())
}

/// Optimizes a list of chunks into a list of non-overlapping chunks.
///
/// Unlike [`merge_chunks`], `min_offset` (typically from the linear index) is given to remove
/// chunks that cannot be in the query.
///
/// # Examples
///
/// ```
/// use noodles_bgzf as bgzf;
/// use noodles_csi::{
///     binning_index::optimize_chunks,
///     index::reference_sequence::bin::Chunk,
/// };
///
/// let chunks = [
///     Chunk::new(bgzf::VirtualPosition::from(2), bgzf::VirtualPosition::from(3)),
///     Chunk::new(bgzf::VirtualPosition::from(5), bgzf::VirtualPosition::from(8)),
///     Chunk::new(bgzf::VirtualPosition::from(7), bgzf::VirtualPosition::from(13)),
///     Chunk::new(bgzf::VirtualPosition::from(21), bgzf::VirtualPosition::from(34)),
/// ];
/// let min_offset = bgzf::VirtualPosition::from(5);
///
/// let actual = optimize_chunks(&chunks, min_offset);
///
/// let expected = [
///     Chunk::new(bgzf::VirtualPosition::from(5), bgzf::VirtualPosition::from(13)),
///     Chunk::new(bgzf::VirtualPosition::from(21), bgzf::VirtualPosition::from(34)),
/// ];
///
/// assert_eq!(actual, expected);
/// ```
pub fn optimize_chunks(chunks: &[Chunk], min_offset: bgzf::VirtualPosition) -> Vec<Chunk> {
    let mut chunks: Vec<_> = chunks
        .iter()
        .filter(|c| c.end() > min_offset)
        .copied()
        .collect();

    if chunks.is_empty() {
        return chunks;
    }

    chunks.sort_unstable_by_key(|c| c.start());

    // At worst, no chunks are merged, and the resulting list will be the same size as the input.
    let mut merged_chunks = Vec::with_capacity(chunks.len());

    // `chunks` is guaranteed to be non-empty.
    let mut current_chunk = chunks[0];

    for next_chunk in chunks.iter().skip(1) {
        if next_chunk.start() > current_chunk.end() {
            merged_chunks.push(current_chunk);
            current_chunk = *next_chunk;
        } else if current_chunk.end() < next_chunk.end() {
            current_chunk = Chunk::new(current_chunk.start(), next_chunk.end());
        }
    }

    merged_chunks.push(current_chunk);

    merged_chunks
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_chunks() -> Vec<Chunk> {
        vec![
            Chunk::new(
                bgzf::VirtualPosition::from(2),
                bgzf::VirtualPosition::from(5),
            ),
            Chunk::new(
                bgzf::VirtualPosition::from(3),
                bgzf::VirtualPosition::from(4),
            ),
            Chunk::new(
                bgzf::VirtualPosition::from(5),
                bgzf::VirtualPosition::from(7),
            ),
            Chunk::new(
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(12),
            ),
            Chunk::new(
                bgzf::VirtualPosition::from(10),
                bgzf::VirtualPosition::from(15),
            ),
            Chunk::new(
                bgzf::VirtualPosition::from(16),
                bgzf::VirtualPosition::from(21),
            ),
        ]
    }

    #[test]
    fn test_merge_chunks() {
        let chunks = build_chunks();
        let actual = merge_chunks(&chunks);

        let expected = [
            Chunk::new(
                bgzf::VirtualPosition::from(2),
                bgzf::VirtualPosition::from(7),
            ),
            Chunk::new(
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(15),
            ),
            Chunk::new(
                bgzf::VirtualPosition::from(16),
                bgzf::VirtualPosition::from(21),
            ),
        ];

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_merge_chunks_with_empty_list() {
        let chunks = Vec::new();
        let merged_chunks = merge_chunks(&chunks);
        assert!(merged_chunks.is_empty());
    }

    #[test]
    fn test_optimize_chunks() {
        let chunks = build_chunks();
        let actual = optimize_chunks(&chunks, bgzf::VirtualPosition::from(10));

        let expected = [
            Chunk::new(
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(15),
            ),
            Chunk::new(
                bgzf::VirtualPosition::from(16),
                bgzf::VirtualPosition::from(21),
            ),
        ];

        assert_eq!(actual, expected);
    }
}
