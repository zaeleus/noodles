//! Binning index reference sequence.

pub mod bin;
pub mod index;
mod metadata;

pub use self::{bin::Bin, index::Index, metadata::Metadata};

use std::io;

use bit_vec::BitVec;
use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use noodles_core::{Position, region::Interval};

use self::bin::Chunk;
use super::resolve_interval;
use crate::binning_index;

/// A binning index reference sequence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceSequence<I> {
    bins: IndexMap<usize, Bin>,
    index: I,
    metadata: Option<Metadata>,
}

impl<I> ReferenceSequence<I>
where
    I: Index,
{
    /// Creates a binning index reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::binning_index::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Default::default(), Vec::new(), None);
    /// ```
    pub fn new(bins: IndexMap<usize, Bin>, index: I, metadata: Option<Metadata>) -> Self {
        Self {
            bins,
            index,
            metadata,
        }
    }

    /// Returns the list of bins in the reference sequence.
    ///
    /// This list does not include the metadata pseudo-bin. Use
    /// [`binning_index::ReferenceSequence::metadata`] instead.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::binning_index::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Default::default(), Vec::new(), None);
    /// assert!(reference_sequence.bins().is_empty());
    /// ```
    pub fn bins(&self) -> &IndexMap<usize, Bin> {
        &self.bins
    }

    /// Returns the index.
    ///
    /// The index is optional and can be empty.
    pub fn index(&self) -> &I {
        &self.index
    }

    /// Returns a list of bins in this reference sequence that intersects the given range.
    ///
    /// The interval values are 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_csi::binning_index::index::ReferenceSequence;
    ///
    /// let reference_sequence = ReferenceSequence::new(Default::default(), Vec::new(), None);
    /// let start = Position::try_from(8)?;
    /// let end = Position::try_from(13)?;
    ///
    /// let query_bins = reference_sequence.query(14, 5, start..=end)?;
    /// assert!(query_bins.is_empty());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<J>(&self, min_shift: u8, depth: u8, interval: J) -> io::Result<Vec<&Bin>>
    where
        J: Into<Interval>,
    {
        let (start, end) = resolve_interval(min_shift, depth, interval)?;

        let max_bin_id = Bin::max_id(depth);
        let mut region_bins = BitVec::from_elem(max_bin_id, false);

        reg2bins(start, end, min_shift, depth, &mut region_bins);

        let query_bins = self
            .bins()
            .iter()
            .filter(|(id, _)| region_bins[**id])
            .map(|(_, bin)| bin)
            .collect();

        Ok(query_bins)
    }

    /// Finds the start virtual position of the first record in the bin that contains that given
    /// start position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_core::Position;
    /// use noodles_csi::binning_index::index::{
    ///     reference_sequence::{index::BinnedIndex, Bin},
    ///     ReferenceSequence,
    /// };
    ///
    /// const MIN_SHIFT: u8 = 4;
    /// const DEPTH: u8 = 2;
    ///
    /// let bins = [(1, Bin::new(Vec::new()))].into_iter().collect();
    /// let index = [(1, bgzf::VirtualPosition::from(233))].into_iter().collect();
    /// let reference_sequence: ReferenceSequence<BinnedIndex> = ReferenceSequence::new(bins, index, None);
    ///
    /// let start = Position::try_from(8)?;
    /// assert_eq!(
    ///     reference_sequence.min_offset(MIN_SHIFT, DEPTH, start),
    ///     bgzf::VirtualPosition::from(233)
    /// );
    ///
    /// let start = Position::try_from(144)?;
    /// assert_eq!(
    ///     reference_sequence.min_offset(MIN_SHIFT, DEPTH, start),
    ///     bgzf::VirtualPosition::default()
    /// );
    ///
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn min_offset(&self, min_shift: u8, depth: u8, start: Position) -> bgzf::VirtualPosition {
        self.index.min_offset(min_shift, depth, start)
    }

    /// Returns the start position of the first record in the last linear bin.
    ///
    /// This uses the linear index, if available; otherwise, the largest linear offset of all the
    /// bins is returned. If there is neither a linear index nor a collection of bins, this returns
    /// `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::binning_index::index::{
    ///     reference_sequence::{index::BinnedIndex, Bin},
    ///     ReferenceSequence,
    /// };
    ///
    /// let reference_sequence = ReferenceSequence::new(
    ///     Default::default(),
    ///     vec![
    ///         bgzf::VirtualPosition::from(8),
    ///         bgzf::VirtualPosition::from(13),
    ///         bgzf::VirtualPosition::from(21),
    ///     ],
    ///     None,
    /// );
    /// assert_eq!(
    ///     reference_sequence.first_record_in_last_linear_bin_start_position(),
    ///     Some(bgzf::VirtualPosition::from(21))
    /// );
    ///
    /// let mut reference_sequence: ReferenceSequence<BinnedIndex> = ReferenceSequence::new(
    ///     [
    ///         (0, Bin::new(Vec::new())),
    ///         (2, Bin::new(Vec::new())),
    ///         (9, Bin::new(Vec::new())),
    ///     ]
    ///     .into_iter()
    ///     .collect(),
    ///     [
    ///         (0, bgzf::VirtualPosition::from(8)),
    ///         (2, bgzf::VirtualPosition::from(21)),
    ///         (9, bgzf::VirtualPosition::from(13)),
    ///     ]
    ///     .into_iter()
    ///     .collect(),
    ///     None,
    /// );
    ///
    /// assert_eq!(
    ///     reference_sequence.first_record_in_last_linear_bin_start_position(),
    ///     Some(bgzf::VirtualPosition::from(21))
    /// );
    ///
    /// let reference_sequence = ReferenceSequence::new(Default::default(), Vec::new(), None);
    /// assert!(reference_sequence.first_record_in_last_linear_bin_start_position().is_none());
    /// ```
    pub fn first_record_in_last_linear_bin_start_position(&self) -> Option<bgzf::VirtualPosition> {
        self.index.last_first_start_position()
    }

    pub(crate) fn update(
        &mut self,
        min_shift: u8,
        depth: u8,
        start: Position,
        end: Position,
        is_mapped: bool,
        chunk: Chunk,
    ) {
        let id = reg2bin(start, end, min_shift, depth);
        let bins = self.bins.entry(id).or_insert(Bin::new(Vec::new()));
        bins.add_chunk(chunk);

        self.index.update(min_shift, depth, start, end, chunk);

        let metadata = self.metadata.get_or_insert(Metadata::new(
            bgzf::VirtualPosition::MAX,
            bgzf::VirtualPosition::MIN,
            0,
            0,
        ));

        metadata.update(is_mapped, chunk);
    }

    pub(crate) fn optimize(mut self) -> Self {
        const MIN_COMPRESSED_SIZE: u64 = 1 << 16;

        let mut bins = IndexMap::with_capacity(self.bins.len());

        let mut ids: Vec<_> = self.bins.keys().copied().collect();
        ids.sort_unstable();

        for id in ids.iter().rev() {
            // SAFETY: The set of IDs is equivalent to the set of bin keys.
            let bin = self.bins.swap_remove(id).unwrap();

            let Some(pid) = parent_id(*id) else {
                bins.insert(*id, bin);
                continue;
            };

            let size = bin_compressed_span(&bin);

            if size < MIN_COMPRESSED_SIZE
                && let Some(parent_bin) = self.bins.get_mut(&pid)
            {
                parent_bin.chunks.extend(bin.chunks);
            } else {
                bins.insert(*id, bin);
            }
        }

        for bin in bins.values_mut() {
            bin.chunks = merge_chunks(bin.chunks());
        }

        Self::new(bins, self.index, self.metadata)
    }
}

fn bin_compressed_span(bin: &Bin) -> u64 {
    let chunks = bin.chunks();

    let (first_chunk, last_chunk) = if chunks.is_empty() {
        return 0;
    } else if chunks.len() == 1 {
        let chunk = chunks[0];
        (chunk, chunk)
    } else {
        chunks_minmax(chunks)
    };

    let start = first_chunk.start().compressed();
    let end = last_chunk.end().compressed();

    end - start
}

fn chunks_minmax(chunks: &[Chunk]) -> (Chunk, Chunk) {
    let chunk = &chunks[0];
    let (mut min, mut max) = (chunk, chunk);

    for chunk in chunks.iter().skip(1) {
        if chunk.start() < min.start() {
            min = chunk;
        }

        if chunk.start() > max.start() {
            max = chunk;
        }
    }

    (*min, *max)
}

fn merge_chunks(chunks: &[Chunk]) -> Vec<Chunk> {
    let mut chunks = chunks.to_vec();

    if chunks.is_empty() {
        return chunks;
    }

    chunks.sort_unstable_by_key(|c| c.start());

    // At worst, no chunks are merged, and the resulting list will be the same size as the input.
    let mut merged_chunks = Vec::with_capacity(chunks.len());

    // `chunks` is guaranteed to be non-empty.
    let mut current_chunk = chunks[0];

    for next_chunk in chunks.iter().skip(1) {
        if current_chunk.end().compressed() >= next_chunk.start().compressed() {
            if current_chunk.end() < next_chunk.end() {
                current_chunk = Chunk::new(current_chunk.start(), next_chunk.end());
            }
        } else {
            merged_chunks.push(current_chunk);
            current_chunk = *next_chunk;
        }
    }

    merged_chunks.push(current_chunk);

    merged_chunks
}

impl<I> binning_index::ReferenceSequence for ReferenceSequence<I>
where
    I: Index,
{
    fn metadata(&self) -> Option<&Metadata> {
        self.metadata.as_ref()
    }
}

const M: usize = 8;

// parent of i = floor((i - 1) / M)
fn parent_id(id: usize) -> Option<usize> {
    if id > 0 { Some((id - 1) / M) } else { None }
}

// `CSIv1.pdf` (2020-07-21)
fn reg2bin(start: Position, end: Position, min_shift: u8, depth: u8) -> usize {
    // [beg, end), 0-based
    let beg = usize::from(start) - 1;
    let end = usize::from(end);

    let end = end - 1;
    let mut l = depth;
    let mut s = min_shift;
    let mut t = ((1 << (depth * 3)) - 1) / 7;

    while l > 0 {
        if beg >> s == end >> s {
            return t + (beg >> s);
        }

        l -= 1;
        s += 3;
        t -= 1 << (l * 3);
    }

    0
}

// `CSIv1.pdf` (2020-07-21)
#[allow(clippy::many_single_char_names)]
fn reg2bins(start: Position, end: Position, min_shift: u8, depth: u8, bins: &mut BitVec) {
    // [beg, end), 0-based
    let beg = usize::from(start) - 1;
    let end = usize::from(end);

    let end = end - 1;
    let mut l = 0;
    let mut t = 0;
    let mut s = i32::from(min_shift) + i32::from(depth) * 3;

    while l <= depth {
        let b = t + (beg >> s);
        let e = t + (end >> s);

        for i in b..=e {
            bins.set(i, true);
        }

        s -= 3;
        t += 1 << (l * 3);
        l += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MIN_SHIFT: u8 = 14;
    const DEPTH: u8 = 5;

    #[cfg(not(target_pointer_width = "16"))]
    #[test]
    fn test_query() {
        let reference_sequence = ReferenceSequence::new(Default::default(), Vec::new(), None);
        let end = const { Position::new(i32::MAX as usize).unwrap() };

        assert!(matches!(
            reference_sequence.query(MIN_SHIFT, DEPTH, ..=end),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));
    }

    #[test]
    fn test_optimize() {
        let bins = [
            (0, vec![((0, 0), (16384, 0))]),
            (1, vec![((0, 0), (32768, 0))]),
            (2, vec![((0, 0), (131072, 0))]),
            (9, vec![((0, 0), (16384, 0))]),
            (10, vec![((0, 0), (131072, 0))]),
        ]
        .into_iter()
        .map(|(id, raw_chunks)| {
            let chunks = raw_chunks
                .into_iter()
                .map(|((a, b), (c, d))| {
                    Chunk::new(
                        bgzf::VirtualPosition::new(a, b).unwrap(),
                        bgzf::VirtualPosition::new(c, d).unwrap(),
                    )
                })
                .collect();

            (id, Bin::new(chunks))
        })
        .collect();

        let reference_sequence = ReferenceSequence::new(bins, Vec::new(), None);
        let actual = reference_sequence.optimize();

        let expected_bins: IndexMap<_, _> = [
            (
                10,
                Bin::new(vec![Chunk::new(
                    const { bgzf::VirtualPosition::new(0, 0).unwrap() },
                    const { bgzf::VirtualPosition::new(131072, 0).unwrap() },
                )]),
            ),
            (
                2,
                Bin::new(vec![Chunk::new(
                    const { bgzf::VirtualPosition::new(0, 0).unwrap() },
                    const { bgzf::VirtualPosition::new(131072, 0).unwrap() },
                )]),
            ),
            (
                0,
                Bin::new(vec![Chunk::new(
                    const { bgzf::VirtualPosition::new(0, 0).unwrap() },
                    const { bgzf::VirtualPosition::new(32768, 0).unwrap() },
                )]),
            ),
        ]
        .into_iter()
        .collect();

        assert_eq!(actual.bins, expected_bins);
    }

    #[test]
    fn test_chunks_minmax() {
        let min_chunk = Chunk::new(
            const { bgzf::VirtualPosition::new(5, 0).unwrap() },
            const { bgzf::VirtualPosition::new(8, 0).unwrap() },
        );

        let max_chunk = Chunk::new(
            const { bgzf::VirtualPosition::new(21, 0).unwrap() },
            const { bgzf::VirtualPosition::new(34, 0).unwrap() },
        );

        let chunks = [
            max_chunk,
            min_chunk,
            Chunk::new(
                const { bgzf::VirtualPosition::new(13, 0).unwrap() },
                const { bgzf::VirtualPosition::new(21, 0).unwrap() },
            ),
        ];

        assert_eq!(chunks_minmax(&chunks), (min_chunk, max_chunk));
    }

    #[test]
    fn test_reg2bin() {
        const MIN_SHIFT: u8 = 4;
        const DEPTH: u8 = 2;

        let start = const { Position::new(8).unwrap() };
        let end = start;
        assert_eq!(reg2bin(start, end, MIN_SHIFT, DEPTH), 9);

        let end = const { Position::new(13).unwrap() };
        assert_eq!(reg2bin(start, end, MIN_SHIFT, DEPTH), 9);

        let end = const { Position::new(16).unwrap() };
        assert_eq!(reg2bin(start, end, MIN_SHIFT, DEPTH), 9);

        let end = const { Position::new(17).unwrap() };
        assert_eq!(reg2bin(start, end, MIN_SHIFT, DEPTH), 1);

        let end = const { Position::new(143).unwrap() };
        assert_eq!(reg2bin(start, end, MIN_SHIFT, DEPTH), 0);
    }

    #[test]
    fn test_reg2bins() -> Result<(), noodles_core::position::TryFromIntError> {
        // +------------------------------------------------------------------------------------...
        // | 0                                                                                  ...
        // | 0-1023                                                                             ...
        // +-------------------------------------------------------------------------+----------...
        // | 1                                                                       | 2        ...
        // | 0-127                                                                   | 128-255  ...
        // +--------+--------+--------+--------+--------+--------+---------+---------+---------+...
        // | 9      | 10     | 11     | 12     | 13     | 14     | 15      | 16      | 17      |...
        // | 0-15   | 16-31  | 32-47  | 48-63  | 64-79  | 80-95  | 96-111  | 112-127 | 128-143 |...
        // +--------+--------+--------+--------+--------+--------+---------+---------+---------+...

        const MIN_SHIFT: u8 = 4;
        const DEPTH: u8 = 2;

        fn t(start: Position, end: Position, expected_bin_ids: &[usize]) {
            let max_bin_id = Bin::max_id(DEPTH);

            let mut actual = BitVec::from_elem(max_bin_id, false);
            reg2bins(start, end, MIN_SHIFT, DEPTH, &mut actual);

            let mut expected = BitVec::from_elem(max_bin_id, false);

            for &bin_id in expected_bin_ids {
                expected.set(bin_id, true);
            }

            assert_eq!(actual, expected);
        }

        t(
            const { Position::new(1).unwrap() },
            const { Position::new(16).unwrap() },
            &[0, 1, 9],
        );

        t(
            const { Position::new(9).unwrap() },
            const { Position::new(13).unwrap() },
            &[0, 1, 9],
        );

        t(
            const { Position::new(36).unwrap() },
            const { Position::new(67).unwrap() },
            &[0, 1, 11, 12, 13],
        );

        t(
            const { Position::new(49).unwrap() },
            const { Position::new(143).unwrap() },
            &[0, 1, 2, 12, 13, 14, 15, 16, 17],
        );

        Ok(())
    }
}
