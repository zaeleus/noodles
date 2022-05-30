mod builder;

pub use self::builder::Builder;

use std::io;

use noodles_core::Position;

use crate::{
    header::{ReferenceSequence, ReferenceSequences},
    record::{Cigar, Data, Flags, MappingQuality, QualityScores, ReadName, Sequence},
    AlignmentRecord,
};

/// An alignment record.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    read_name: Option<ReadName>,
    flags: Flags,
    reference_sequence_id: Option<usize>,
    alignment_start: Option<Position>,
    mapping_quality: Option<MappingQuality>,
    cigar: Cigar,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<Position>,
    template_length: i32,
    sequence: Sequence,
    quality_scores: QualityScores,
    data: Data,
}

impl Record {
    /// Creates an alignment record builder.
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns a mutable reference to the read name.
    pub fn read_name_mut(&mut self) -> &mut Option<ReadName> {
        &mut self.read_name
    }

    /// Returns a mutable reference to the flags.
    pub fn flags_mut(&mut self) -> &mut Flags {
        &mut self.flags
    }

    /// Returns the reference sequence ID.
    pub fn reference_sequence_id(&self) -> Option<usize> {
        self.reference_sequence_id
    }

    /// Returns a mutable reference to the reference sequence ID.
    pub fn reference_sequence_id_mut(&mut self) -> &mut Option<usize> {
        &mut self.reference_sequence_id
    }

    /// Returns a mutable reference to the alignment start.
    pub fn alignment_start_mut(&mut self) -> &mut Option<Position> {
        &mut self.alignment_start
    }

    /// Returns a mutable reference to the mapping quality.
    pub fn mapping_quality_mut(&mut self) -> &mut Option<MappingQuality> {
        &mut self.mapping_quality
    }

    /// Returns a mutable reference to the CIGAR operations.
    pub fn cigar_mut(&mut self) -> &mut Cigar {
        &mut self.cigar
    }

    /// Returns the mate reference sequence ID.
    pub fn mate_reference_sequence_id(&self) -> Option<usize> {
        self.mate_reference_sequence_id
    }

    /// Returns a mutable reference to the mate reference sequence ID.
    pub fn mate_reference_sequence_id_mut(&mut self) -> &mut Option<usize> {
        &mut self.mate_reference_sequence_id
    }

    /// Returns a mutable reference to the mate alignment start.
    pub fn mate_alignment_start_mut(&mut self) -> &mut Option<Position> {
        &mut self.mate_alignment_start
    }

    /// Returns a mutable reference to template length.
    pub fn template_length_mut(&mut self) -> &mut i32 {
        &mut self.template_length
    }

    /// Returns a mutable reference to sequence.
    pub fn sequence_mut(&mut self) -> &mut Sequence {
        &mut self.sequence
    }

    /// Returns a mutable reference to quality scores.
    pub fn quality_scores_mut(&mut self) -> &mut QualityScores {
        &mut self.quality_scores
    }

    /// Returns a mutable reference to the data.
    pub fn data_mut(&mut self) -> &mut Data {
        &mut self.data
    }
}

impl AlignmentRecord for Record {
    fn read_name(&self) -> Option<&ReadName> {
        self.read_name.as_ref()
    }

    fn reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<std::io::Result<&'rs ReferenceSequence>> {
        get_reference_sequence(reference_sequences, self.reference_sequence_id())
    }

    fn flags(&self) -> Flags {
        self.flags
    }

    fn alignment_start(&self) -> Option<Position> {
        self.alignment_start
    }

    fn alignment_span(&self) -> usize {
        self.cigar().reference_len()
    }

    fn mapping_quality(&self) -> Option<MappingQuality> {
        self.mapping_quality
    }

    fn cigar(&self) -> &Cigar {
        &self.cigar
    }

    fn mate_reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<std::io::Result<&'rs ReferenceSequence>> {
        get_reference_sequence(reference_sequences, self.mate_reference_sequence_id())
    }

    fn mate_alignment_start(&self) -> Option<Position> {
        self.mate_alignment_start
    }

    fn template_length(&self) -> i32 {
        self.template_length
    }

    fn sequence(&self) -> &Sequence {
        &self.sequence
    }

    fn quality_scores(&self) -> &QualityScores {
        &self.quality_scores
    }

    fn data(&self) -> &Data {
        &self.data
    }
}

fn get_reference_sequence(
    reference_sequences: &ReferenceSequences,
    reference_sequence_id: Option<usize>,
) -> Option<io::Result<&ReferenceSequence>> {
    reference_sequence_id.map(|id| {
        reference_sequences
            .get_index(id)
            .map(|(_, rs)| rs)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "invalid reference sequence ID")
            })
    })
}
