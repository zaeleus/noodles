use std::{cmp, collections::HashMap};

use noodles_bgzf as bgzf;

use crate::Record;

use super::{
    bin::{self, Chunk},
    Bin, Metadata, ReferenceSequence,
};

const WINDOW_SIZE: i32 = 16384;

// § 5.1.1 Basic binning index: "... bins 4681-37448 span 16Kbp regions."
const MAX_INTERVAL_COUNT: usize = 37448 - 4681 + 1;

#[derive(Debug)]
pub struct Builder {
    bin_builders: HashMap<u32, bin::Builder>,
    intervals: Vec<bgzf::VirtualPosition>,
    start_position: bgzf::VirtualPosition,
    end_position: bgzf::VirtualPosition,
    mapped_record_count: u64,
    unmapped_record_count: u64,
}

impl Builder {
    pub fn add_record(&mut self, record: &Record, chunk: Chunk) {
        self.update_metadata(record, chunk);

        let bin_id = record.bin() as u32;

        let builder = self.bin_builders.entry(bin_id).or_insert_with(|| {
            let mut builder = Bin::builder();
            builder.set_id(bin_id);
            builder
        });

        builder.add_chunk(chunk);

        let start = i32::from(record.position());
        let reference_len = record.cigar().reference_len() as i32;
        let end = start + reference_len - 1;

        let linear_index_start_offset = ((start - 1) / WINDOW_SIZE) as usize;
        let linear_index_end_offset = ((end - 1) / WINDOW_SIZE) as usize;

        for i in linear_index_start_offset..=linear_index_end_offset {
            self.intervals[i] = chunk.start();
        }
    }

    pub fn build(self) -> ReferenceSequence {
        let bins: Vec<_> = self
            .bin_builders
            .into_iter()
            .map(|(_, b)| b.build())
            .collect();

        let metadata = Metadata::new(
            self.start_position,
            self.end_position,
            self.mapped_record_count,
            self.unmapped_record_count,
        );

        ReferenceSequence::new(bins, self.intervals, Some(metadata))
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
            intervals: vec![bgzf::VirtualPosition::default(); MAX_INTERVAL_COUNT],
            start_position: bgzf::VirtualPosition::max(),
            end_position: bgzf::VirtualPosition::default(),
            mapped_record_count: 0,
            unmapped_record_count: 0,
        }
    }
}
