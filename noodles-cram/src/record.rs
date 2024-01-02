//! CRAM record and fields.

mod builder;
mod convert;
pub mod feature;
mod features;
mod flags;
mod next_mate_flags;
pub mod resolve;

pub use self::{
    builder::Builder, feature::Feature, features::Features, flags::Flags,
    next_mate_flags::NextMateFlags,
};

use std::io;

use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::record_buf::{QualityScores, Sequence},
    header::record::value::{
        map::{self, ReferenceSequence},
        Map,
    },
};

/// A CRAM record.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    pub(crate) id: u64,
    pub(crate) bam_bit_flags: sam::record::Flags,
    pub(crate) cram_bit_flags: Flags,
    pub(crate) reference_sequence_id: Option<usize>,
    pub(crate) read_length: usize,
    pub(crate) alignment_start: Option<Position>,
    pub(crate) read_group_id: Option<usize>,
    pub(crate) name: Option<sam::record::Name>,
    pub(crate) next_mate_bit_flags: NextMateFlags,
    pub(crate) next_fragment_reference_sequence_id: Option<usize>,
    pub(crate) next_mate_alignment_start: Option<Position>,
    pub(crate) template_size: i32,
    pub(crate) distance_to_next_fragment: Option<usize>,
    pub(crate) tags: sam::record::Data,
    pub(crate) bases: Sequence,
    pub(crate) features: Features,
    pub(crate) mapping_quality: Option<sam::record::MappingQuality>,
    pub(crate) quality_scores: QualityScores,
}

impl Record {
    /// Returns a builder to create a record from each of its fields.
    pub fn builder() -> Builder {
        Builder::default()
    }

    pub(crate) fn id(&self) -> u64 {
        self.id
    }

    /// Returns the BAM flags.
    ///
    /// This is also called the BAM bit flags.
    pub fn bam_flags(&self) -> sam::record::Flags {
        self.bam_bit_flags
    }

    /// Returns the SAM flags.
    pub fn flags(&self) -> sam::record::Flags {
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
    pub fn reference_sequence_id(&self) -> Option<usize> {
        self.reference_sequence_id
    }

    /// Returns the associated reference sequence.
    pub fn reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs sam::header::ReferenceSequences,
    ) -> Option<
        io::Result<(
            &'rs map::reference_sequence::Name,
            &'rs Map<ReferenceSequence>,
        )>,
    > {
        get_reference_sequence(reference_sequences, self.reference_sequence_id())
    }

    /// Returns the read length.
    pub fn read_length(&self) -> usize {
        self.read_length
    }

    /// Returns the alignment start.
    pub fn alignment_start(&self) -> Option<Position> {
        self.alignment_start
    }

    /// Returns the alignment span.
    fn alignment_span(&self) -> usize {
        calculate_alignment_span(self.read_length(), self.features())
    }

    /// Returns the alignment end.
    pub fn alignment_end(&self) -> Option<Position> {
        self.alignment_start().and_then(|alignment_start| {
            let end = usize::from(alignment_start) + self.alignment_span() - 1;
            Position::new(end)
        })
    }

    /// Returns the read group ID.
    ///
    /// This is also simply called the read group. It is the position of the read group in the SAM
    /// header.
    pub fn read_group_id(&self) -> Option<usize> {
        self.read_group_id
    }

    /// Returns the name.
    pub fn name(&self) -> Option<&sam::record::Name> {
        self.name.as_ref()
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
    pub fn next_fragment_reference_sequence_id(&self) -> Option<usize> {
        self.next_fragment_reference_sequence_id
    }

    /// Returns the associated mate reference sequence.
    pub fn mate_reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs sam::header::ReferenceSequences,
    ) -> Option<
        io::Result<(
            &'rs map::reference_sequence::Name,
            &'rs Map<ReferenceSequence>,
        )>,
    > {
        get_reference_sequence(
            reference_sequences,
            self.next_fragment_reference_sequence_id(),
        )
    }

    /// Returns the alignment start position of the next mate.
    ///
    /// This value is 1-based.
    pub fn next_mate_alignment_start(&self) -> Option<Position> {
        self.next_mate_alignment_start
    }

    /// Returns the alignment start.
    pub fn mate_alignment_start(&self) -> Option<Position> {
        self.next_mate_alignment_start
    }

    /// Returns the template size.
    pub fn template_size(&self) -> i32 {
        self.template_size
    }

    /// Returns the template size.
    pub fn template_length(&self) -> i32 {
        self.template_size
    }

    /// Returns the distance to the next fragment.
    ///
    /// This is the number of records to the next fragment within a slice.
    pub fn distance_to_next_fragment(&self) -> Option<usize> {
        self.distance_to_next_fragment
    }

    /// Returns the tag dictionary.
    pub fn tags(&self) -> &sam::record::Data {
        &self.tags
    }

    /// Returns the data.
    pub fn data(&self) -> &sam::record::Data {
        &self.tags
    }

    /// Returns the read bases.
    pub fn bases(&self) -> &Sequence {
        &self.bases
    }

    /// Returns the sequence.
    pub fn sequence(&self) -> &Sequence {
        &self.bases
    }

    /// Returns the read features.
    pub fn features(&self) -> &Features {
        &self.features
    }

    /// Returns the mapping quality.
    pub fn mapping_quality(&self) -> Option<sam::record::MappingQuality> {
        self.mapping_quality
    }

    /// Returns the quality scores.
    pub fn quality_scores(&self) -> &QualityScores {
        &self.quality_scores
    }
}

impl Default for Record {
    fn default() -> Self {
        Builder::default().build()
    }
}

pub(crate) fn calculate_alignment_span(read_length: usize, features: &Features) -> usize {
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
    reference_sequence_id: Option<usize>,
) -> Option<io::Result<(&map::reference_sequence::Name, &Map<ReferenceSequence>)>> {
    reference_sequence_id.map(|id| {
        reference_sequences.get_index(id).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "invalid reference sequence ID")
        })
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_alignment_span() -> Result<(), noodles_core::position::TryFromIntError> {
        let features = Features::default();
        assert_eq!(calculate_alignment_span(4, &features), 4);

        let features = Features::from(vec![Feature::HardClip(Position::try_from(1)?, 4)]);
        assert_eq!(calculate_alignment_span(4, &features), 4);

        let features = Features::from(vec![
            Feature::Insertion(Position::try_from(1)?, vec![b'A', b'C']),
            Feature::InsertBase(Position::try_from(4)?, b'G'),
            Feature::Deletion(Position::try_from(6)?, 3),
            Feature::ReferenceSkip(Position::try_from(10)?, 5),
            Feature::SoftClip(Position::try_from(16)?, vec![b'A', b'C', b'G', b'T']),
        ]);
        assert_eq!(calculate_alignment_span(20, &features), 21);

        Ok(())
    }
}
