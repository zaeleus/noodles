use std::io;

use indexmap::IndexMap;
use noodles_core::Position;

use super::index::{
    reference_sequence::{self, bin::Chunk},
    Header, Index, ReferenceSequence,
};

/// A binning index indexer.
#[derive(Debug)]
pub struct Indexer<I> {
    min_shift: u8,
    depth: u8,
    header: Option<Header>,
    reference_sequences: Vec<ReferenceSequence<I>>,
    unplaced_unmapped_record_count: u64,
}

impl<I> Indexer<I>
where
    I: reference_sequence::Index + Default,
{
    /// Creates a binning index indexer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::binning_index::{index::reference_sequence::index::BinnedIndex, Indexer};
    /// let indexer = Indexer::<BinnedIndex>::new(14, 5);
    /// ```
    pub fn new(min_shift: u8, depth: u8) -> Self {
        Self {
            min_shift,
            depth,
            header: None,
            reference_sequences: Vec::new(),
            unplaced_unmapped_record_count: 0,
        }
    }

    /// Sets a tabix header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::binning_index::{
    ///     index::{reference_sequence::index::BinnedIndex, Header},
    ///     Indexer,
    /// };
    ///
    /// let header = Header::default();
    /// let indexer = Indexer::<BinnedIndex>::new(14, 5).set_header(header);
    /// ```
    pub fn set_header(mut self, header: Header) -> Self {
        self.header = Some(header);
        self
    }

    /// Adds a record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_core::Position;
    /// use noodles_csi::binning_index::{
    ///     index::reference_sequence::{bin::Chunk, index::BinnedIndex},
    ///     Indexer,
    /// };
    ///
    /// let mut indexer = Indexer::<BinnedIndex>::new(14, 5);
    ///
    /// let reference_sequence_id = 0;
    /// let start = Position::try_from(8)?;
    /// let end = Position::try_from(13)?;
    /// let is_mapped = true;
    /// let chunk = Chunk::new(
    ///     bgzf::VirtualPosition::from(144),
    ///     bgzf::VirtualPosition::from(233),
    /// );
    ///
    /// indexer.add_record(Some((reference_sequence_id, start, end, is_mapped)), chunk)?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn add_record(
        &mut self,
        alignment_context: Option<(usize, Position, Position, bool)>,
        chunk: Chunk,
    ) -> io::Result<()> {
        use std::cmp::Ordering;

        let Some((reference_sequence_id, start, end, is_mapped)) = alignment_context else {
            self.unplaced_unmapped_record_count += 1;
            return Ok(());
        };

        if self.reference_sequences.is_empty() {
            self.add_reference_sequences_until(0);
        }

        // SAFETY: `self.reference_sequences` is non-empty.
        let current_reference_sequence_id = self.reference_sequences.len() - 1;

        match reference_sequence_id.cmp(&current_reference_sequence_id) {
            Ordering::Less => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "invalid reference sequence ID",
                ));
            }
            Ordering::Equal => {}
            Ordering::Greater => self.add_reference_sequences_until(reference_sequence_id),
        }

        let reference_sequence = &mut self.reference_sequences[reference_sequence_id];
        reference_sequence.update(self.min_shift, self.depth, start, end, is_mapped, chunk);

        Ok(())
    }

    /// Builds a binning index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::binning_index::{index::reference_sequence::index::BinnedIndex, Indexer};
    /// let indexer = Indexer::<BinnedIndex>::new(14, 5);
    /// let index = indexer.build(0);
    /// ```
    pub fn build(mut self, reference_sequence_count: usize) -> Index<I> {
        if reference_sequence_count == 0 {
            return Index::builder()
                .set_unplaced_unmapped_record_count(self.unplaced_unmapped_record_count)
                .build();
        }

        // SAFETY: `reference_sequence_count` is > 0.
        self.add_reference_sequences_until(reference_sequence_count - 1);

        let mut builder = Index::builder()
            .set_reference_sequences(self.reference_sequences)
            .set_unplaced_unmapped_record_count(self.unplaced_unmapped_record_count);

        if let Some(header) = self.header {
            builder = builder.set_header(header);
        }

        builder.build()
    }

    fn add_reference_sequences_until(&mut self, reference_sequence_id: usize) {
        let new_len = reference_sequence_id + 1;

        self.reference_sequences.resize_with(new_len, || {
            ReferenceSequence::new(IndexMap::new(), Default::default(), None)
        });
    }
}

impl<I> Default for Indexer<I>
where
    I: reference_sequence::Index + Default,
{
    fn default() -> Self {
        Self {
            min_shift: 14,
            depth: 5,
            header: None,
            reference_sequences: Vec::new(),
            unplaced_unmapped_record_count: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use noodles_bgzf as bgzf;

    use super::*;
    use crate::binning_index::index::reference_sequence::{index::LinearIndex, Bin, Metadata};

    #[test]
    fn test_default() {
        let indexer = Indexer::<LinearIndex>::default();

        assert_eq!(indexer.min_shift, 14);
        assert_eq!(indexer.depth, 5);
        assert!(indexer.header.is_none());
        assert!(indexer.reference_sequences.is_empty());
        assert_eq!(indexer.unplaced_unmapped_record_count, 0);
    }

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        const MIN_SHIFT: u8 = 14;
        const DEPTH: u8 = 5;

        let mut indexer = Indexer::<LinearIndex>::new(MIN_SHIFT, DEPTH);

        indexer.add_record(
            Some((0, Position::try_from(8)?, Position::try_from(13)?, true)),
            Chunk::new(
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(9),
            ),
        )?;

        indexer.add_record(
            Some((
                0,
                Position::try_from(121393)?,
                Position::try_from(196418)?,
                false,
            )),
            Chunk::new(
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(3473408),
            ),
        )?;

        let actual = indexer.build(1);

        let expected = {
            let bins = [
                (
                    4681,
                    Bin::new(vec![Chunk::new(
                        bgzf::VirtualPosition::from(0),
                        bgzf::VirtualPosition::from(9),
                    )]),
                ),
                (
                    73,
                    Bin::new(vec![Chunk::new(
                        bgzf::VirtualPosition::from(9),
                        bgzf::VirtualPosition::from(3473408),
                    )]),
                ),
            ]
            .into_iter()
            .collect();

            let index = vec![
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
            ];

            let metadata = Metadata::new(
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(3473408),
                1,
                1,
            );

            Index::builder()
                .set_reference_sequences(vec![ReferenceSequence::new(bins, index, Some(metadata))])
                .set_unplaced_unmapped_record_count(0)
                .build()
        };

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_build_with_reference_sequence_count() {
        let index = Indexer::<LinearIndex>::default().build(2);
        assert_eq!(index.reference_sequences().len(), 2);
    }
}
