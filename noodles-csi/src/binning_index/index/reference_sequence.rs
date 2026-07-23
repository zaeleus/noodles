//! Binning index reference sequence.

pub mod bin;
pub mod index;
mod metadata;

pub use self::{bin::Bin, index::Index, metadata::Metadata};

use std::{io, num::NonZero};

use bit_vec::BitVec;
use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use noodles_core::{Position, region::Interval};

use self::bin::Chunk;
use super::resolve_interval;
use crate::binning_index;

// A bin is kept only if its records span at least this many bytes of compressed data; a bin
// covering less is folded into its parent during compaction. Same value as htslib's
// HTS_MIN_MARKER_DIST (`compress_binning`, hts.c).
const MIN_COMPRESSED_SPAN: u64 = 1 << 16;

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

    // Folds bins that span little compressed data into their parents and merges chunks that share a
    // BGZF block, matching htslib's `compress_binning`. This significantly reduces the size of the
    // index without changing query results.
    pub(crate) fn compact(&mut self) {
        // Snapshot the insertion order so serialization stays deterministic after the folding
        // below (bins are an ordered map for this reason).
        let order: Vec<usize> = self.bins.keys().copied().collect();

        // First pass: fold a bin into its parent when its records span less than
        // `MIN_COMPRESSED_SPAN` bytes of compressed data. Bins are visited deepest first
        // (descending ID) so that a bin has received all of its descendants' folded chunks before
        // its own span is tested, and so a parent is never removed before its children are
        // considered.
        let mut ids = order.clone();
        ids.sort_unstable();

        for &id in ids.iter().rev() {
            let Some(pid) = parent_id(id) else { continue };

            let Some(bin) = self.bins.get_mut(&id) else {
                continue;
            };

            bin.sort_chunks();

            if bin.compressed_span() >= MIN_COMPRESSED_SPAN {
                continue;
            }

            if !self.bins.contains_key(&pid) {
                continue;
            }

            if let Some(bin) = self.bins.swap_remove(&id)
                && let Some(parent) = self.bins.get_mut(&pid)
            {
                parent.extend_chunks(bin.into_chunks());
            }
        }

        // Second pass: restore the original insertion order and, within each surviving bin, sort by
        // start and merge chunks that share a BGZF block.
        let mut compacted = IndexMap::with_capacity(self.bins.len());

        for id in order {
            if let Some(mut bin) = self.bins.swap_remove(&id) {
                bin.sort_chunks();
                bin.merge_chunks_in_same_block();
                compacted.insert(id, bin);
            }
        }

        self.bins = compacted;
    }
}

impl<I> binning_index::ReferenceSequence for ReferenceSequence<I>
where
    I: Index,
{
    fn metadata(&self) -> Option<&Metadata> {
        self.metadata.as_ref()
    }
}

const M: usize = NonZero::new(8).unwrap().get();

// parent of i = floor((i - 1) / M)
pub(crate) fn parent_id(id: usize) -> Option<usize> {
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
    use super::index::LinearIndex;
    use super::*;

    const MIN_SHIFT: u8 = 14;
    const DEPTH: u8 = 5;

    // Builds a virtual position from a compressed (block) offset and an uncompressed offset.
    fn vp(compressed: u64, uncompressed: u16) -> bgzf::VirtualPosition {
        bgzf::VirtualPosition::try_from((compressed, uncompressed)).unwrap()
    }

    #[cfg(not(target_pointer_width = "16"))]
    #[test]
    fn test_query() -> Result<(), noodles_core::position::TryFromIntError> {
        let reference_sequence = ReferenceSequence::new(Default::default(), Vec::new(), None);
        let end = usize::try_from(i32::MAX).and_then(Position::try_from)?;

        assert!(matches!(
            reference_sequence.query(MIN_SHIFT, DEPTH, ..=end),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }

    #[test]
    fn test_reg2bin() -> Result<(), noodles_core::position::TryFromIntError> {
        const MIN_SHIFT: u8 = 4;
        const DEPTH: u8 = 2;

        let start = Position::try_from(8)?;
        let end = start;
        assert_eq!(reg2bin(start, end, MIN_SHIFT, DEPTH), 9);

        let end = Position::try_from(13)?;
        assert_eq!(reg2bin(start, end, MIN_SHIFT, DEPTH), 9);

        let end = Position::try_from(16)?;
        assert_eq!(reg2bin(start, end, MIN_SHIFT, DEPTH), 9);

        let end = Position::try_from(17)?;
        assert_eq!(reg2bin(start, end, MIN_SHIFT, DEPTH), 1);

        let end = Position::try_from(143)?;
        assert_eq!(reg2bin(start, end, MIN_SHIFT, DEPTH), 0);

        Ok(())
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

        t(Position::try_from(1)?, Position::try_from(16)?, &[0, 1, 9]);
        t(Position::try_from(9)?, Position::try_from(13)?, &[0, 1, 9]);

        t(
            Position::try_from(36)?,
            Position::try_from(67)?,
            &[0, 1, 11, 12, 13],
        );

        t(
            Position::try_from(49)?,
            Position::try_from(143)?,
            &[0, 1, 2, 12, 13, 14, 15, 16, 17],
        );

        Ok(())
    }

    #[test]
    fn test_compact_folds_bin_into_parent() {
        // Bin 4681 is a child of bin 585 (parent_id(4681) == 585).
        let bins = [
            (585, Bin::new(vec![Chunk::new(vp(1, 0), vp(1, 100))])),
            (4681, Bin::new(vec![Chunk::new(vp(1, 100), vp(1, 200))])),
        ]
        .into_iter()
        .collect();
        let mut reference_sequence: ReferenceSequence<LinearIndex> =
            ReferenceSequence::new(bins, Vec::new(), None);

        reference_sequence.compact();

        assert_eq!(reference_sequence.bins().len(), 1);
        assert_eq!(
            reference_sequence.bins()[&585].chunks(),
            [Chunk::new(vp(1, 0), vp(1, 200))]
        );
    }

    #[test]
    fn test_compact_keeps_bin_without_parent() {
        let bins = [(4681, Bin::new(vec![Chunk::new(vp(1, 0), vp(1, 100))]))]
            .into_iter()
            .collect();
        let mut reference_sequence: ReferenceSequence<LinearIndex> =
            ReferenceSequence::new(bins, Vec::new(), None);

        reference_sequence.compact();

        assert!(reference_sequence.bins().contains_key(&4681));
    }

    #[test]
    fn test_compact_keeps_bin_with_large_span() {
        let bins = [
            (585, Bin::new(vec![Chunk::new(vp(0, 0), vp(0, 10))])),
            (4681, Bin::new(vec![Chunk::new(vp(0, 0), vp(1 << 16, 0))])),
        ]
        .into_iter()
        .collect();
        let mut reference_sequence: ReferenceSequence<LinearIndex> =
            ReferenceSequence::new(bins, Vec::new(), None);

        reference_sequence.compact();

        assert!(reference_sequence.bins().contains_key(&585));
        assert!(reference_sequence.bins().contains_key(&4681));
    }

    #[test]
    fn test_compact_does_not_fold_root() {
        let bins = [(0, Bin::new(vec![Chunk::new(vp(0, 0), vp(0, 100))]))]
            .into_iter()
            .collect();
        let mut reference_sequence: ReferenceSequence<LinearIndex> =
            ReferenceSequence::new(bins, Vec::new(), None);

        reference_sequence.compact();

        assert!(reference_sequence.bins().contains_key(&0));
    }

    #[test]
    fn test_compact_cascades() {
        // Bin 4681 -> 585 -> 73 (each is the parent of the next); bin 9 (73's parent) is absent.
        let bins = [
            (73, Bin::new(vec![Chunk::new(vp(1, 0), vp(1, 10))])),
            (585, Bin::new(vec![Chunk::new(vp(1, 10), vp(1, 20))])),
            (4681, Bin::new(vec![Chunk::new(vp(1, 20), vp(1, 30))])),
        ]
        .into_iter()
        .collect();
        let mut reference_sequence: ReferenceSequence<LinearIndex> =
            ReferenceSequence::new(bins, Vec::new(), None);

        reference_sequence.compact();

        assert_eq!(reference_sequence.bins().len(), 1);
        assert_eq!(
            reference_sequence.bins()[&73].chunks(),
            [Chunk::new(vp(1, 0), vp(1, 30))]
        );
    }
}
