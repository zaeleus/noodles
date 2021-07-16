use std::{cmp, collections::HashMap};

use noodles_bgzf as bgzf;
use noodles_csi::index::reference_sequence::bin::Chunk;

use super::{bin, Bin, Metadata, ReferenceSequence, WINDOW_SIZE};

#[derive(Debug, Default)]
pub struct Builder {
    bin_builders: HashMap<u32, bin::Builder>,
    intervals: Vec<Option<bgzf::VirtualPosition>>,
    start_position: bgzf::VirtualPosition,
    end_position: bgzf::VirtualPosition,
    mapped_record_count: u64,
}

impl Builder {
    pub fn add_record(&mut self, start: i32, end: i32, chunk: Chunk) -> &mut Self {
        self.update_bins(start, end, chunk);
        self.update_linear_index(start, end, chunk);
        self.update_metadata(chunk);
        self
    }

    pub fn build(self) -> ReferenceSequence {
        if self.bin_builders.is_empty() {
            return ReferenceSequence::default();
        }

        let bins = self
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
            0,
        );

        ReferenceSequence::new(bins, intervals, Some(metadata))
    }

    fn update_bins(&mut self, start: i32, end: i32, chunk: Chunk) {
        let bin_id = region_to_bin(start, end) as u32;

        let builder = self.bin_builders.entry(bin_id).or_insert_with(|| {
            let mut builder = Bin::builder();
            builder.set_id(bin_id);
            builder
        });

        builder.add_chunk(chunk);
    }

    fn update_linear_index(&mut self, start: i32, end: i32, chunk: Chunk) {
        let linear_index_start_offset = ((start - 1) / WINDOW_SIZE) as usize;
        let linear_index_end_offset = ((end - 1) / WINDOW_SIZE) as usize;

        if linear_index_end_offset >= self.intervals.len() {
            self.intervals
                .resize(linear_index_end_offset + 1, Default::default());
        }

        for i in linear_index_start_offset..=linear_index_end_offset {
            self.intervals[i].get_or_insert(chunk.start());
        }
    }

    fn update_metadata(&mut self, chunk: Chunk) {
        self.mapped_record_count += 1;
        self.start_position = cmp::min(self.start_position, chunk.start());
        self.end_position = cmp::max(self.end_position, chunk.end());
    }
}

// 0-based, [start, end)
#[allow(clippy::eq_op)]
fn region_to_bin(start: i32, mut end: i32) -> i32 {
    end -= 1;

    if start >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (start >> 14)
    } else if start >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (start >> 17)
    } else if start >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (start >> 20)
    } else if start >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (start >> 23)
    } else if start >> 26 == end >> 26 {
        ((1 << 3) - 1) / 7 + (start >> 26)
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build() {
        let mut builder = Builder::default();

        builder.add_record(
            8,
            13,
            Chunk::new(
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(9),
            ),
        );

        builder.add_record(
            121393,
            196418,
            Chunk::new(
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(3473408),
            ),
        );

        let actual = builder.build();

        let expected = {
            let bins = vec![
                Bin::new(
                    4681,
                    vec![Chunk::new(
                        bgzf::VirtualPosition::from(0),
                        bgzf::VirtualPosition::from(9),
                    )],
                ),
                Bin::new(
                    73,
                    vec![Chunk::new(
                        bgzf::VirtualPosition::from(9),
                        bgzf::VirtualPosition::from(3473408),
                    )],
                ),
            ];

            let intervals = vec![
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
            ];

            ReferenceSequence::new(bins, intervals, None)
        };

        for expected_bin in expected.bins() {
            let actual_bin = actual
                .bins()
                .iter()
                .find(|b| b.id() == expected_bin.id())
                .expect("missing bin");

            assert_eq!(actual_bin, expected_bin);
        }

        assert_eq!(actual.intervals(), expected.intervals());
    }

    #[test]
    fn test_build_with_no_bins() {
        let reference_sequence = Builder::default().build();
        assert_eq!(reference_sequence, ReferenceSequence::default());
    }

    #[test]
    fn test_region_to_bin() {
        // [8, 13]
        assert_eq!(region_to_bin(7, 13), 4681);
        // [63245986, 63245986]
        assert_eq!(region_to_bin(63245985, 63255986), 8541);
    }
}
