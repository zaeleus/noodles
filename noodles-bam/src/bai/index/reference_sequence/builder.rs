use std::{cmp, collections::HashMap, io};

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::index::reference_sequence::bin::Chunk;
use noodles_sam::record::Flags;

use crate::writer::record::region_to_bin;

use super::{bin, Bin, Metadata, ReferenceSequence, MIN_SHIFT};

// § 5.2 The BAI index format for BAM files (2020-07-19)
const MAX_INTERVAL_COUNT: usize = 131072;

#[derive(Debug)]
pub struct Builder {
    bin_builders: HashMap<usize, bin::Builder>,
    intervals: Vec<Option<bgzf::VirtualPosition>>,
    start_position: bgzf::VirtualPosition,
    end_position: bgzf::VirtualPosition,
    mapped_record_count: u64,
    unmapped_record_count: u64,
}

impl Builder {
    pub fn add_record(
        &mut self,
        start: Position,
        end: Position,
        flags: Flags,
        chunk: Chunk,
    ) -> io::Result<()> {
        self.update_bins(start, end, chunk)?;
        self.update_linear_index(start, end, chunk);
        self.update_metadata(flags, chunk);
        Ok(())
    }

    pub fn build(self) -> ReferenceSequence {
        if self.bin_builders.is_empty() {
            return ReferenceSequence::default();
        }

        let bins: Vec<_> = self
            .bin_builders
            .into_iter()
            .map(|(_, b)| b.build())
            .collect();

        let intervals = self
            .intervals
            .into_iter()
            .map(|p| p.unwrap_or_default())
            .collect();

        let metadata = Metadata::new(
            self.start_position,
            self.end_position,
            self.mapped_record_count,
            self.unmapped_record_count,
        );

        ReferenceSequence::new(bins, intervals, Some(metadata))
    }

    fn update_bins(&mut self, start: Position, end: Position, chunk: Chunk) -> io::Result<()> {
        let bin_id = region_to_bin(start, end).map(usize::from)?;

        let builder = self.bin_builders.entry(bin_id).or_insert_with(|| {
            let mut builder = Bin::builder();
            builder.set_id(bin_id);
            builder
        });

        builder.add_chunk(chunk);

        Ok(())
    }

    fn update_linear_index(&mut self, start: Position, end: Position, chunk: Chunk) {
        const WINDOW_SIZE: usize = 1 << MIN_SHIFT;

        let linear_index_start_offset = (usize::from(start) - 1) / WINDOW_SIZE;
        let linear_index_end_offset = (usize::from(end) - 1) / WINDOW_SIZE;

        if linear_index_end_offset >= self.intervals.len() {
            self.intervals
                .resize(linear_index_end_offset + 1, Default::default());
        }

        for i in linear_index_start_offset..=linear_index_end_offset {
            self.intervals[i].get_or_insert(chunk.start());
        }
    }

    fn update_metadata(&mut self, flags: Flags, chunk: Chunk) {
        if flags.is_unmapped() {
            self.unmapped_record_count += 1;
        } else {
            self.mapped_record_count += 1;
        }

        self.start_position = cmp::min(self.start_position, chunk.start());
        self.end_position = cmp::max(self.end_position, chunk.end());
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            bin_builders: HashMap::new(),
            intervals: Vec::with_capacity(MAX_INTERVAL_COUNT),
            start_position: bgzf::VirtualPosition::max(),
            end_position: bgzf::VirtualPosition::default(),
            mapped_record_count: 0,
            unmapped_record_count: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        let mut builder = Builder::default();

        builder.add_record(
            Position::try_from(2)?,
            Position::try_from(5)?,
            Flags::empty(),
            Chunk::new(
                bgzf::VirtualPosition::from(55),
                bgzf::VirtualPosition::from(89),
            ),
        )?;

        builder.add_record(
            Position::try_from(6)?,
            Position::try_from(7)?,
            Flags::UNMAPPED,
            Chunk::new(
                bgzf::VirtualPosition::from(89),
                bgzf::VirtualPosition::from(144),
            ),
        )?;

        let actual = builder.build();

        let expected = ReferenceSequence::new(
            vec![Bin::new(
                4681,
                vec![Chunk::new(
                    bgzf::VirtualPosition::from(55),
                    bgzf::VirtualPosition::from(144),
                )],
            )],
            vec![bgzf::VirtualPosition::from(55)],
            Some(Metadata::new(
                bgzf::VirtualPosition::from(55),
                bgzf::VirtualPosition::from(144),
                1,
                1,
            )),
        );

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_build_with_no_bins() {
        let reference_sequence = Builder::default().build();
        assert_eq!(reference_sequence, ReferenceSequence::default());
    }

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert!(builder.bin_builders.is_empty());
        assert!(builder.intervals.is_empty());

        assert_eq!(builder.start_position, bgzf::VirtualPosition::max());
        assert_eq!(builder.end_position, bgzf::VirtualPosition::default());
        assert_eq!(builder.mapped_record_count, 0);
        assert_eq!(builder.unmapped_record_count, 0);
    }
}
