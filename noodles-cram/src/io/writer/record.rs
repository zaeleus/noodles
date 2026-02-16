mod convert;
mod feature;

use bstr::BString;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        record::{MappingQuality, data::field::Tag},
        record_buf::data::field::Value,
    },
};

pub use self::feature::Feature;
use crate::record::{Flags, MateFlags};

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Record {
    pub(crate) bam_flags: sam::alignment::record::Flags,
    pub(crate) cram_flags: Flags,
    pub(crate) reference_sequence_id: Option<usize>,
    pub(crate) read_length: usize,
    pub(crate) alignment_start: Option<Position>,
    pub(crate) read_group_id: Option<usize>,
    pub(crate) name: Option<BString>,
    pub(crate) mate_flags: MateFlags,
    pub(crate) mate_reference_sequence_id: Option<usize>,
    pub(crate) mate_alignment_start: Option<Position>,
    pub(crate) template_length: i64,
    pub(crate) mate_distance: Option<usize>,
    pub(crate) data: Vec<(Tag, Value)>,
    pub(crate) features: Vec<Feature>,
    pub(crate) mapping_quality: Option<MappingQuality>,
    pub(crate) sequence: Vec<u8>,
    pub(crate) quality_scores: Vec<u8>,
}

impl Record {
    pub fn alignment_span(&self) -> usize {
        calculate_alignment_span(self.read_length, &self.features)
    }

    pub fn alignment_end(&self) -> Option<Position> {
        self.alignment_start.and_then(|start| {
            let end = usize::from(start) + self.alignment_span() - 1;
            Position::new(end)
        })
    }
}

fn calculate_alignment_span(read_length: usize, features: &[Feature]) -> usize {
    features
        .iter()
        .fold(read_length, |span, feature| match feature {
            Feature::Insertion { bases, .. } => span - bases.len(),
            Feature::InsertBase { .. } => span - 1,
            Feature::Deletion { len, .. } => span + len,
            Feature::ReferenceSkip { len, .. } => span + len,
            Feature::SoftClip { bases, .. } => span - bases.len(),
            _ => span,
        })
}
