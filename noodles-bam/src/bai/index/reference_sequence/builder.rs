use std::collections::HashMap;

use noodles_bgzf as bgzf;

use crate::Record;

use super::{
    bin::{self, Chunk},
    Bin, ReferenceSequence,
};

const WINDOW_SIZE: i32 = 16384;

#[derive(Debug)]
pub struct Builder {
    bin_builders: HashMap<u32, bin::Builder>,
    intervals: Vec<bgzf::VirtualPosition>,
}

impl Builder {
    pub fn add_record(&mut self, record: &Record, chunk: Chunk) {
        let bin_id = record.bin() as u32;

        let builder = self.bin_builders.entry(bin_id).or_insert_with(|| {
            let mut builder = Bin::builder();
            builder.set_bin(bin_id);
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
        let bins = self
            .bin_builders
            .into_iter()
            .map(|(_, b)| b.build())
            .collect();

        ReferenceSequence {
            bins,
            intervals: self.intervals,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            bin_builders: HashMap::new(),
            intervals: vec![bgzf::VirtualPosition::default(); 40000],
        }
    }
}
