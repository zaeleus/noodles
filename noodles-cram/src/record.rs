pub mod feature;
mod flags;
mod next_mate_flags;
pub mod resolve;
pub mod tag;

pub use self::{feature::Feature, flags::Flags, next_mate_flags::NextMateFlags, tag::Tag};

use std::{fmt, str};

use noodles_sam as sam;

#[derive(Clone, Default)]
pub struct Record {
    pub id: i64,
    pub bam_bit_flags: sam::record::Flags,
    pub cram_bit_flags: i32,
    pub reference_id: i32,
    pub read_length: i32,
    pub alignment_start: i32,
    pub read_group: i32,
    pub read_name: Vec<u8>,
    pub next_mate_bit_flags: i32,
    pub next_fragment_reference_sequence_id: i32,
    pub next_mate_alignment_start: i32,
    pub template_size: i32,
    pub distance_to_next_fragment: i32,
    pub tags: Vec<Tag>,
    pub bases: Vec<u8>,
    pub features: Vec<Feature>,
    pub mapping_quality: i32,
    pub quality_scores: Vec<u8>,
}

impl Record {
    pub fn bam_bit_flags(&self) -> sam::record::Flags {
        self.bam_bit_flags
    }

    pub fn cram_bit_flags(&self) -> Flags {
        // `cram_bit_flags` can safely be casted to to a u8 because CRAM currently only has 4 bit
        // flags, i.e., the largest value is 2^4 - 1.
        Flags::from(self.cram_bit_flags as u8)
    }

    pub fn read_length(&self) -> i32 {
        self.read_length
    }

    pub fn alignment_start(&self) -> i32 {
        self.alignment_start
    }

    pub fn next_mate_bit_flags(&self) -> NextMateFlags {
        // `next_mate_bit_flags` can safely be casted to a u8 because the `MF` data series only has
        // 2 bit flags, i.e., the largest value is 2^2 - 1;
        NextMateFlags::from(self.next_mate_bit_flags as u8)
    }

    pub fn add_tag(&mut self, tag: Tag) {
        self.tags.push(tag);
    }

    pub fn features(&self) -> &[Feature] {
        &self.features
    }

    pub fn add_feature(&mut self, feature: Feature) {
        self.features.push(feature);
    }
}

impl fmt::Debug for Record {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        let read_name = str::from_utf8(&self.read_name);

        fmt.debug_struct("Record")
            .field("id", &self.id)
            .field("bam_bit_flags", &self.bam_bit_flags())
            .field("cram_bit_flags", &self.cram_bit_flags())
            .field("reference_id", &self.reference_id)
            .field("read_length", &self.read_length)
            .field("alignment_start", &self.alignment_start)
            .field("read_group", &self.read_group)
            .field("read_name", &read_name)
            .field("next_mate_bit_flags", &self.next_mate_bit_flags())
            .field(
                "next_fragment_reference_sequence_id",
                &self.next_fragment_reference_sequence_id,
            )
            .field("next_mate_alignment_start", &self.next_mate_alignment_start)
            .field("template_size", &self.template_size)
            .field("distance_to_next_fragment", &self.distance_to_next_fragment)
            .field("tags", &self.tags)
            .field("bases", &self.bases)
            .field("features", &self.features)
            .field("mapping_quality", &self.mapping_quality)
            .field("quality_scores", &self.quality_scores)
            .finish()
    }
}
