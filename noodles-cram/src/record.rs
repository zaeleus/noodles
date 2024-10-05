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

use bstr::{BStr, BString};
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        record::MappingQuality,
        record_buf::{Data, QualityScores, Sequence},
    },
    header::record::value::{map::ReferenceSequence, Map},
};

/// A CRAM record.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    pub(crate) id: u64,
    pub(crate) bam_bit_flags: sam::alignment::record::Flags,
    pub(crate) cram_bit_flags: Flags,
    pub(crate) reference_sequence_id: Option<usize>,
    pub(crate) read_length: usize,
    pub(crate) alignment_start: Option<Position>,
    pub(crate) read_group_id: Option<usize>,
    pub(crate) name: Option<BString>,
    pub(crate) next_mate_bit_flags: NextMateFlags,
    pub(crate) next_fragment_reference_sequence_id: Option<usize>,
    pub(crate) next_mate_alignment_start: Option<Position>,
    pub(crate) template_size: i32,
    pub(crate) distance_to_next_fragment: Option<usize>,
    pub(crate) tags: Data,
    pub(crate) bases: Sequence,
    pub(crate) features: Features,
    pub(crate) mapping_quality: Option<MappingQuality>,
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
    pub fn bam_flags(&self) -> sam::alignment::record::Flags {
        self.bam_bit_flags
    }

    /// Returns the SAM flags.
    pub fn flags(&self) -> sam::alignment::record::Flags {
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
    pub fn reference_sequence<'h>(
        &self,
        reference_sequences: &'h sam::header::ReferenceSequences,
    ) -> Option<io::Result<(&'h [u8], &'h Map<ReferenceSequence>)>> {
        get_reference_sequence(reference_sequences, self.reference_sequence_id())
    }

    /// Returns the read length.
    pub fn read_length(&self) -> usize {
        self.read_length
    }

    /// Returns the alignment start.
    ///
    /// This position is 1-based, inclusive.
    pub fn alignment_start(&self) -> Option<Position> {
        self.alignment_start
    }

    /// Returns the alignment span.
    fn alignment_span(&self) -> usize {
        calculate_alignment_span(self.read_length(), self.features())
    }

    /// Returns the alignment end.
    ///
    /// This position is 1-based, inclusive.
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
    pub fn name(&self) -> Option<&BStr> {
        self.name.as_ref().map(|name| name.as_ref())
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
    pub fn mate_reference_sequence<'h>(
        &self,
        reference_sequences: &'h sam::header::ReferenceSequences,
    ) -> Option<io::Result<(&'h [u8], &'h Map<ReferenceSequence>)>> {
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
    pub fn tags(&self) -> &Data {
        &self.tags
    }

    /// Returns the data.
    pub fn data(&self) -> &Data {
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
    pub fn mapping_quality(&self) -> Option<MappingQuality> {
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

struct Cigar<'a> {
    features: &'a Features,
    is_unmapped: bool,
    read_length: usize,
}

impl<'a> Cigar<'a> {
    fn new(features: &'a Features, is_unmapped: bool, read_length: usize) -> Self {
        Self {
            features,
            is_unmapped,
            read_length,
        }
    }
}

impl<'a> sam::alignment::record::Cigar for Cigar<'a> {
    fn is_empty(&self) -> bool {
        self.is_unmapped
    }

    fn len(&self) -> usize {
        if self.is_unmapped {
            0
        } else {
            self.iter().count()
        }
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<sam::alignment::record::cigar::Op>> + '_> {
        use std::iter;

        use sam::alignment::record::cigar::iter::TrySimplify;

        if self.is_unmapped {
            Box::new(iter::empty())
        } else {
            Box::new(TrySimplify::new(
                self.features.cigar(self.read_length).map(Ok),
            ))
        }
    }
}

impl sam::alignment::Record for Record {
    fn name(&self) -> Option<&BStr> {
        self.name()
    }

    fn flags(&self) -> io::Result<sam::alignment::record::Flags> {
        Ok(self.flags())
    }

    fn reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        _: &'h sam::Header,
    ) -> Option<io::Result<usize>> {
        self.reference_sequence_id().map(Ok)
    }

    fn alignment_start(&self) -> Option<io::Result<Position>> {
        self.alignment_start().map(Ok)
    }

    fn mapping_quality(&self) -> Option<io::Result<MappingQuality>> {
        self.mapping_quality().map(Ok)
    }

    fn cigar(&self) -> Box<dyn sam::alignment::record::Cigar + '_> {
        Box::new(Cigar::new(
            &self.features,
            self.flags().is_unmapped(),
            self.read_length(),
        ))
    }

    fn mate_reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        _: &'h sam::Header,
    ) -> Option<io::Result<usize>> {
        self.next_fragment_reference_sequence_id().map(Ok)
    }

    fn mate_alignment_start(&self) -> Option<io::Result<Position>> {
        self.next_mate_alignment_start().map(Ok)
    }

    fn template_length(&self) -> io::Result<i32> {
        Ok(self.template_length())
    }

    fn sequence(&self) -> Box<dyn sam::alignment::record::Sequence + '_> {
        Box::new(self.sequence())
    }

    fn quality_scores(&self) -> Box<dyn sam::alignment::record::QualityScores + '_> {
        Box::new(self.quality_scores())
    }

    fn data(&self) -> Box<dyn sam::alignment::record::Data + '_> {
        Box::new(self.data())
    }
}

pub(crate) fn calculate_alignment_span(read_length: usize, features: &Features) -> usize {
    features
        .iter()
        .fold(read_length, |alignment_span, feature| match feature {
            Feature::Insertion { bases, .. } => alignment_span - bases.len(),
            Feature::InsertBase { .. } => alignment_span - 1,
            Feature::Deletion { len, .. } => alignment_span + len,
            Feature::ReferenceSkip { len, .. } => alignment_span + len,
            Feature::SoftClip { bases, .. } => alignment_span - bases.len(),
            _ => alignment_span,
        })
}

fn get_reference_sequence(
    reference_sequences: &sam::header::ReferenceSequences,
    reference_sequence_id: Option<usize>,
) -> Option<io::Result<(&[u8], &Map<ReferenceSequence>)>> {
    reference_sequence_id.map(|id| {
        reference_sequences
            .get_index(id)
            .map(|(name, map)| (name.as_ref(), map))
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
        let features = Features::default();
        assert_eq!(calculate_alignment_span(4, &features), 4);

        let features = Features::from(vec![Feature::HardClip {
            position: Position::try_from(1)?,
            len: 4,
        }]);
        assert_eq!(calculate_alignment_span(4, &features), 4);

        let features = Features::from(vec![
            Feature::Insertion {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::InsertBase {
                position: Position::try_from(4)?,
                base: b'G',
            },
            Feature::Deletion {
                position: Position::try_from(6)?,
                len: 3,
            },
            Feature::ReferenceSkip {
                position: Position::try_from(10)?,
                len: 5,
            },
            Feature::SoftClip {
                position: Position::try_from(16)?,
                bases: vec![b'A', b'C', b'G', b'T'],
            },
        ]);
        assert_eq!(calculate_alignment_span(20, &features), 21);

        Ok(())
    }
}
