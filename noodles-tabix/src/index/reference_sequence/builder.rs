use std::collections::HashMap;

use super::{
    bin::{self, Chunk},
    Bin, ReferenceSequence,
};

#[derive(Debug, Default)]
pub struct Builder {
    bin_builders: HashMap<u32, bin::Builder>,
}

impl Builder {
    pub fn add_record(&mut self, start: u32, end: u32, chunk: Chunk) -> &mut Self {
        self.update_bins(start, end, chunk);
        self
    }

    pub fn build(self) -> ReferenceSequence {
        let bins = self
            .bin_builders
            .into_iter()
            .map(|(_, b)| b.build())
            .collect();

        // TODO: linear index
        ReferenceSequence::new(bins, Vec::new())
    }

    fn update_bins(&mut self, start: u32, end: u32, chunk: Chunk) {
        let bin_id = region_to_bin(start, end);

        let builder = self.bin_builders.entry(bin_id).or_insert_with(|| {
            let mut builder = Bin::builder();
            builder.set_id(bin_id);
            builder
        });

        builder.add_chunk(chunk);
    }
}

// 0-based, [start, end)
#[allow(clippy::eq_op)]
fn region_to_bin(start: u32, mut end: u32) -> u32 {
    assert!(end > 0);

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
