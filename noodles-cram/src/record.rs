//! CRAM record and fields.

mod data;
pub mod feature;
mod features;
mod flags;
mod mate_flags;
mod quality_scores;
pub mod resolve;

pub use self::{feature::Feature, flags::Flags, mate_flags::MateFlags};

use std::io;

use bstr::{BStr, BString};
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        record::{data::field::Tag, MappingQuality},
        record_buf::{data::field::Value, Sequence},
    },
    header::record::value::{map::ReferenceSequence, Map},
};

use self::{data::Data, quality_scores::QualityScores};

/// A CRAM record.
#[derive(Clone, Debug, PartialEq)]
pub struct Record<'c> {
    pub(crate) id: u64,
    pub(crate) header: Option<&'c sam::Header>,
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
    pub(crate) template_length: i32,
    pub(crate) distance_to_mate: Option<usize>,
    pub(crate) data: Vec<(Tag, Value)>,
    pub(crate) sequence: Sequence,
    pub(crate) features: Vec<Feature>,
    pub(crate) mapping_quality: Option<MappingQuality>,
    pub(crate) quality_scores: &'c [u8],
}

impl Record<'_> {
    pub(crate) fn id(&self) -> u64 {
        self.id
    }

    /// Returns the BAM flags.
    ///
    /// This is also called the BAM bit flags.
    pub fn bam_flags(&self) -> sam::alignment::record::Flags {
        self.bam_flags
    }

    /// Returns the SAM flags.
    pub fn flags(&self) -> sam::alignment::record::Flags {
        self.bam_flags
    }

    /// Returns the CRAM flags.
    ///
    /// This is also called the CRAM bit flags or compression bit flags.
    pub fn cram_flags(&self) -> Flags {
        self.cram_flags
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
    pub fn next_mate_flags(&self) -> MateFlags {
        self.mate_flags
    }

    /// Returns the reference sequence ID of the next fragment.
    ///
    /// It is the position of the reference sequence in the SAM header.
    pub fn next_fragment_reference_sequence_id(&self) -> Option<usize> {
        self.mate_reference_sequence_id
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
        self.mate_alignment_start
    }

    /// Returns the alignment start.
    pub fn mate_alignment_start(&self) -> Option<Position> {
        self.mate_alignment_start
    }

    /// Returns the template size.
    pub fn template_size(&self) -> i32 {
        self.template_length
    }

    /// Returns the template size.
    pub fn template_length(&self) -> i32 {
        self.template_length
    }

    /// Returns the distance to the next fragment.
    ///
    /// This is the number of records to the next fragment within a slice.
    pub fn distance_to_next_fragment(&self) -> Option<usize> {
        self.distance_to_mate
    }

    /// Returns the tag dictionary.
    pub fn tags(&self) -> &[(Tag, Value)] {
        &self.data
    }

    /// Returns the data.
    pub fn data(&self) -> &[(Tag, Value)] {
        &self.data
    }

    /// Returns the read bases.
    #[deprecated(since = "0.70.0", note = "Use `Record::sequence` instead.")]
    pub fn bases(&self) -> &Sequence {
        &self.sequence
    }

    /// Returns the sequence.
    pub fn sequence(&self) -> &Sequence {
        &self.sequence
    }

    /// Returns the read features.
    pub fn features(&self) -> &[Feature] {
        &self.features
    }

    /// Returns the mapping quality.
    pub fn mapping_quality(&self) -> Option<MappingQuality> {
        self.mapping_quality
    }

    /// Returns the quality scores.
    pub fn quality_scores(&self) -> &[u8] {
        self.quality_scores
    }
}

impl Default for Record<'_> {
    fn default() -> Self {
        Self {
            id: 0,
            header: None,
            bam_flags: sam::alignment::record::Flags::UNMAPPED,
            cram_flags: Flags::default(),
            reference_sequence_id: None,
            read_length: 0,
            alignment_start: None,
            read_group_id: None,
            name: None,
            mate_flags: MateFlags::default(),
            mate_reference_sequence_id: None,
            mate_alignment_start: None,
            template_length: 0,
            distance_to_mate: None,
            data: Vec::new(),
            sequence: Sequence::default(),
            features: Vec::new(),
            mapping_quality: None,
            quality_scores: &[],
        }
    }
}

struct Cigar<'a> {
    features: &'a [Feature],
    is_unmapped: bool,
    read_length: usize,
}

impl<'a> Cigar<'a> {
    fn new(features: &'a [Feature], is_unmapped: bool, read_length: usize) -> Self {
        Self {
            features,
            is_unmapped,
            read_length,
        }
    }
}

impl sam::alignment::record::Cigar for Cigar<'_> {
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
                features::Cigar::new(self.features, self.read_length).map(Ok),
            ))
        }
    }
}

impl sam::alignment::Record for Record<'_> {
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
        if self.bam_flags.is_unmapped() || self.cram_flags.are_quality_scores_stored_as_array() {
            Box::new(Scores(self.quality_scores))
        } else {
            Box::new(QualityScores::new(&self.features, self.read_length))
        }
    }

    fn data(&self) -> Box<dyn sam::alignment::record::Data + '_> {
        Box::new(Data::new(
            self.header.unwrap(),
            &self.data,
            self.read_group_id,
        ))
    }
}

struct Scores<'c>(&'c [u8]);

impl sam::alignment::record::QualityScores for Scores<'_> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u8>> + '_> {
        Box::new(self.0.iter().copied().map(Ok))
    }
}

pub(crate) fn calculate_alignment_span(read_length: usize, features: &[Feature]) -> usize {
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
        let features = [];
        assert_eq!(calculate_alignment_span(4, &features), 4);

        let features = [Feature::HardClip {
            position: Position::try_from(1)?,
            len: 4,
        }];
        assert_eq!(calculate_alignment_span(4, &features), 4);

        let features = [
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
        ];
        assert_eq!(calculate_alignment_span(20, &features), 21);

        Ok(())
    }
}
