use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use noodles_core::Position;

use super::{bin::Chunk, Bin, Index, Metadata, ReferenceSequence};

/// A CSI reference sequence builder.
#[derive(Debug)]
pub struct Builder<I> {
    bins: IndexMap<usize, Bin>,
    index: I,
    start_position: bgzf::VirtualPosition,
    end_position: bgzf::VirtualPosition,
    mapped_record_count: u64,
    unmapped_record_count: u64,
}

impl<I> Builder<I>
where
    I: Index,
{
    /// Adds a record.
    pub fn add_record(
        &mut self,
        min_shift: u8,
        depth: u8,
        start: Position,
        end: Position,
        is_mapped: bool,
        chunk: Chunk,
    ) {
        self.update_bins(min_shift, depth, start, end, chunk);
        self.index.update(min_shift, depth, start, end, chunk);
        self.update_metadata(is_mapped, chunk);
    }

    /// Builds a CSI reference sequence.
    pub fn build(self) -> ReferenceSequence<I> {
        if self.bins.is_empty() {
            return ReferenceSequence::new(IndexMap::new(), self.index, None);
        }

        let metadata = Metadata::new(
            self.start_position,
            self.end_position,
            self.mapped_record_count,
            self.unmapped_record_count,
        );

        ReferenceSequence::new(self.bins, self.index, Some(metadata))
    }

    fn update_bins(
        &mut self,
        min_shift: u8,
        depth: u8,
        start: Position,
        end: Position,
        chunk: Chunk,
    ) {
        use super::reg2bin;

        let bin_id = reg2bin(start, end, min_shift, depth);
        let builder = self.bins.entry(bin_id).or_insert(Bin::new(Vec::new()));
        builder.add_chunk(chunk);
    }

    fn update_metadata(&mut self, is_mapped: bool, chunk: Chunk) {
        if is_mapped {
            self.mapped_record_count += 1;
        } else {
            self.unmapped_record_count += 1;
        }

        self.start_position = self.start_position.min(chunk.start());
        self.end_position = self.end_position.max(chunk.end());
    }
}

impl<I> Default for Builder<I>
where
    I: Index + Default,
{
    fn default() -> Self {
        Self {
            bins: IndexMap::new(),
            index: I::default(),
            start_position: bgzf::VirtualPosition::MAX,
            end_position: bgzf::VirtualPosition::MIN,
            mapped_record_count: 0,
            unmapped_record_count: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binning_index::{
        index::reference_sequence::index::LinearIndex, ReferenceSequence as _,
    };

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        const MIN_SHIFT: u8 = 14;
        const DEPTH: u8 = 5;

        let mut builder = Builder::<LinearIndex>::default();

        builder.add_record(
            MIN_SHIFT,
            DEPTH,
            Position::try_from(8)?,
            Position::try_from(13)?,
            true,
            Chunk::new(
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(9),
            ),
        );

        builder.add_record(
            MIN_SHIFT,
            DEPTH,
            Position::try_from(121393)?,
            Position::try_from(196418)?,
            false,
            Chunk::new(
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(3473408),
            ),
        );

        let actual = builder.build();

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

            ReferenceSequence::new(bins, index, Some(metadata))
        };

        assert_eq!(actual.bins(), expected.bins());
        assert_eq!(actual.index(), expected.index());
        assert_eq!(actual.metadata(), expected.metadata());

        Ok(())
    }
}
