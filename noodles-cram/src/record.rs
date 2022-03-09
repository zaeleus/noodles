//! CRAM record and fields.

mod builder;
mod convert;
pub mod feature;
mod features;
mod flags;
mod next_mate_flags;
pub(crate) mod read_group_id;
pub mod resolve;
pub mod tag;

pub use self::{
    builder::Builder, feature::Feature, features::Features, flags::Flags,
    next_mate_flags::NextMateFlags, read_group_id::ReadGroupId, tag::Tag,
};

use std::{fmt, io};

use noodles_bam as bam;
use noodles_sam as sam;

/// A CRAM record.
#[derive(Clone, PartialEq)]
pub struct Record {
    pub(crate) id: i64,
    pub(crate) bam_bit_flags: sam::record::Flags,
    pub(crate) cram_bit_flags: Flags,
    pub(crate) reference_sequence_id: Option<bam::record::ReferenceSequenceId>,
    pub(crate) read_length: usize,
    pub(crate) alignment_start: Option<sam::record::Position>,
    pub(crate) read_group: Option<ReadGroupId>,
    pub(crate) read_name: Option<sam::record::ReadName>,
    pub(crate) next_mate_bit_flags: NextMateFlags,
    pub(crate) next_fragment_reference_sequence_id: Option<bam::record::ReferenceSequenceId>,
    pub(crate) next_mate_alignment_start: Option<sam::record::Position>,
    pub(crate) template_size: i32,
    pub(crate) distance_to_next_fragment: i32,
    pub(crate) tags: Vec<Tag>,
    pub(crate) bases: Vec<u8>,
    pub(crate) features: Features,
    pub(crate) mapping_quality: Option<sam::record::MappingQuality>,
    pub(crate) quality_scores: sam::record::QualityScores,
}

impl Record {
    /// Returns a builder to create a record from each of its fields.
    pub fn builder() -> Builder {
        Builder::default()
    }

    pub(crate) fn id(&self) -> i64 {
        self.id
    }

    /// Returns the BAM flags.
    ///
    /// This is also called the BAM bit flags.
    pub fn bam_flags(&self) -> sam::record::Flags {
        self.bam_bit_flags
    }

    /// Returns the CRAM flags.
    ///
    /// This is also called the CRAM bit flags or compression bit flags.
    pub fn cram_flags(&self) -> Flags {
        self.cram_bit_flags
    }

    /// Returns the reference sequence ID.
    ///
    /// This is also called the reference ID. It is the position of the reference sequence in the
    /// SAM header.
    pub fn reference_sequence_id(&self) -> Option<bam::record::ReferenceSequenceId> {
        self.reference_sequence_id
    }

    /// Returns the read length.
    pub fn read_length(&self) -> usize {
        self.read_length
    }

    /// Returns the read group ID.
    ///
    /// This is also simply called the read group. It is the position of the read group in the SAM
    /// header.
    pub fn read_group_id(&self) -> Option<ReadGroupId> {
        self.read_group
    }

    /// Returns the next mate flags.
    ///
    /// This is also call the next mate bit flags.
    pub fn next_mate_flags(&self) -> NextMateFlags {
        self.next_mate_bit_flags
    }

    /// Returns the reference sequence ID of the next fragment.
    ///
    /// It is the position of the reference sequence in the SAM header.
    pub fn next_fragment_reference_sequence_id(&self) -> Option<bam::record::ReferenceSequenceId> {
        self.next_fragment_reference_sequence_id
    }

    /// Returns the alignment start position of the next mate.
    ///
    /// This value is 1-based.
    pub fn next_mate_alignment_start(&self) -> Option<sam::record::Position> {
        self.next_mate_alignment_start
    }

    /// Returns the template size.
    pub fn template_size(&self) -> i32 {
        self.template_size
    }

    /// Returns the distance to the next fragment.
    ///
    /// This is the number of records to the next fragment within a slice.
    pub fn distance_to_next_fragment(&self) -> i32 {
        self.distance_to_next_fragment
    }

    /// Returns the tag dictionary.
    pub fn tags(&self) -> &[Tag] {
        &self.tags
    }

    /// Returns the read bases.
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }

    /// Returns the read features.
    pub fn features(&self) -> &Features {
        &self.features
    }

    pub(crate) fn add_feature(&mut self, feature: Feature) {
        self.features.push(feature);
    }
}

impl Default for Record {
    fn default() -> Self {
        Builder::default().build()
    }
}

impl fmt::Debug for Record {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.debug_struct("Record")
            .field("id", &self.id)
            .field("bam_bit_flags", &self.bam_flags())
            .field("cram_bit_flags", &self.cram_flags())
            .field("reference_id", &self.reference_sequence_id)
            .field("read_length", &self.read_length)
            .field("alignment_start", &self.alignment_start)
            .field("read_group", &self.read_group)
            .field("read_name", &self.read_name)
            .field("next_mate_bit_flags", &self.next_mate_flags())
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

impl sam::AlignmentRecord for Record {
    fn read_name(&self) -> Option<&sam::record::ReadName> {
        self.read_name.as_ref()
    }

    fn reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs sam::header::ReferenceSequences,
    ) -> Option<io::Result<&'rs sam::header::ReferenceSequence>> {
        get_reference_sequence(reference_sequences, self.reference_sequence_id())
    }

    fn flags(&self) -> sam::record::Flags {
        self.bam_bit_flags
    }

    fn alignment_start(&self) -> Option<sam::record::Position> {
        self.alignment_start
    }

    fn alignment_span(&self) -> io::Result<u32> {
        Ok(calculate_alignment_span(self.read_length(), self.features()) as u32)
    }

    fn mapping_quality(&self) -> Option<sam::record::MappingQuality> {
        self.mapping_quality
    }

    fn mate_reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs sam::header::ReferenceSequences,
    ) -> Option<io::Result<&'rs sam::header::ReferenceSequence>> {
        get_reference_sequence(
            reference_sequences,
            self.next_fragment_reference_sequence_id(),
        )
    }

    fn mate_alignment_start(&self) -> Option<sam::record::Position> {
        self.next_mate_alignment_start
    }

    fn template_length(&self) -> i32 {
        self.template_size
    }

    fn quality_scores(&self) -> &sam::record::QualityScores {
        &self.quality_scores
    }
}

fn calculate_alignment_span(read_length: usize, features: &Features) -> usize {
    features
        .iter()
        .fold(read_length, |alignment_span, feature| match feature {
            Feature::Insertion(_, bases) => alignment_span - bases.len(),
            Feature::InsertBase(_, _) => alignment_span - 1,
            Feature::Deletion(_, len) => alignment_span + len,
            Feature::ReferenceSkip(_, len) => alignment_span + len,
            Feature::SoftClip(_, bases) => alignment_span - bases.len(),
            _ => alignment_span,
        })
}

fn get_reference_sequence(
    reference_sequences: &sam::header::ReferenceSequences,
    reference_sequence_id: Option<bam::record::ReferenceSequenceId>,
) -> Option<io::Result<&sam::header::ReferenceSequence>> {
    reference_sequence_id.map(|id| {
        reference_sequences
            .get_index(usize::from(id))
            .map(|(_, rs)| rs)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "invalid reference sequence ID")
            })
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_alignment_span() -> Result<(), noodles_core::position::TryFromIntError> {
        use noodles_core::Position;
        use sam::record::sequence::Base;

        let features = Features::default();
        assert_eq!(calculate_alignment_span(4, &features), 4);

        let features = Features::from(vec![Feature::HardClip(Position::try_from(1)?, 4)]);
        assert_eq!(calculate_alignment_span(4, &features), 4);

        let features = Features::from(vec![
            Feature::Insertion(Position::try_from(1)?, vec![Base::A, Base::C]),
            Feature::InsertBase(Position::try_from(4)?, Base::G),
            Feature::Deletion(Position::try_from(6)?, 3),
            Feature::ReferenceSkip(Position::try_from(10)?, 5),
            Feature::SoftClip(
                Position::try_from(16)?,
                vec![Base::A, Base::C, Base::G, Base::T],
            ),
        ]);
        assert_eq!(calculate_alignment_span(20, &features), 21);

        Ok(())
    }
}
