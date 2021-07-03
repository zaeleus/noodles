use std::{cmp, collections::HashMap, io};

use noodles_bgzf::{self as bgzf, index::Chunk};

use crate::Record;

use super::{bin, Bin, Metadata, ReferenceSequence, WINDOW_SIZE};

// § 5.2 The BAI index format for BAM files (2020-07-19)
const MAX_INTERVAL_COUNT: usize = 131072;

#[derive(Debug)]
pub struct Builder {
    bin_builders: HashMap<u32, bin::Builder>,
    intervals: Vec<Option<bgzf::VirtualPosition>>,
    start_position: bgzf::VirtualPosition,
    end_position: bgzf::VirtualPosition,
    mapped_record_count: u64,
    unmapped_record_count: u64,
}

impl Builder {
    pub fn add_record(&mut self, record: &Record, chunk: Chunk) -> io::Result<()> {
        self.update_bins(record, chunk);
        self.update_linear_index(record, chunk)?;
        self.update_metadata(record, chunk);
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

    fn update_bins(&mut self, record: &Record, chunk: Chunk) {
        let bin_id = u32::from(record.bin());

        let builder = self.bin_builders.entry(bin_id).or_insert_with(|| {
            let mut builder = Bin::builder();
            builder.set_id(bin_id);
            builder
        });

        builder.add_chunk(chunk);
    }

    fn update_linear_index(&mut self, record: &Record, chunk: Chunk) -> io::Result<()> {
        let start = record.position().map(i32::from).expect("missing position");
        let reference_len = record.cigar().reference_len().map(|len| len as i32)?;
        let end = start + reference_len - 1;

        let linear_index_start_offset = ((start - 1) / WINDOW_SIZE) as usize;
        let linear_index_end_offset = ((end - 1) / WINDOW_SIZE) as usize;

        if linear_index_end_offset >= self.intervals.len() {
            self.intervals
                .resize(linear_index_end_offset + 1, Default::default());
        }

        for i in linear_index_start_offset..=linear_index_end_offset {
            self.intervals[i].get_or_insert(chunk.start());
        }

        Ok(())
    }

    fn update_metadata(&mut self, record: &Record, chunk: Chunk) {
        if record.flags().is_unmapped() {
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
    use std::convert::TryFrom;

    use noodles_sam::{
        self as sam,
        record::{Flags, Position},
    };

    use super::*;

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequences = vec![("sq0", 8)]
            .into_iter()
            .map(|(name, len)| {
                sam::header::ReferenceSequence::new(name, len).map(|rs| (name.into(), rs))
            })
            .collect::<Result<_, _>>()?;

        let mut builder = Builder::default();

        let record = Record::try_from_sam_record(
            &reference_sequences,
            &sam::Record::builder()
                .set_flags(Flags::empty())
                .set_position(Position::try_from(2)?)
                .set_cigar("4M".parse()?)
                .build()?,
        )?;

        builder.add_record(
            &record,
            Chunk::new(
                bgzf::VirtualPosition::from(55),
                bgzf::VirtualPosition::from(89),
            ),
        )?;

        let record = Record::try_from_sam_record(
            &reference_sequences,
            &sam::Record::builder()
                .set_position(Position::try_from(6)?)
                .set_cigar("2M".parse()?)
                .build()?,
        )?;

        builder.add_record(
            &record,
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
