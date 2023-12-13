use super::record::{
    Cigar, Data, Flags, MappingQuality, Position, QualityScores, ReadName, ReferenceSequenceId,
    Sequence, TemplateLength,
};
use crate::Header;

/// An alignment record.
pub trait AnyRecord {
    /// Returns the read name.
    fn read_name(&self) -> Option<&dyn ReadName>;

    /// Returns the flags.
    fn flags(&self) -> &dyn Flags;

    /// Returns the reference sequence ID.
    fn reference_sequence_id(&self, header: &Header) -> Option<&dyn ReferenceSequenceId>;

    /// Returns the alignment start.
    fn alignment_start(&self) -> Option<&dyn Position>;

    /// Returns the mapping quality.
    fn mapping_quality(&self) -> Option<&dyn MappingQuality>;

    /// Returns the CIGAR operations.
    fn cigar(&self, header: &Header) -> &dyn Cigar;

    /// Returns the mate reference sequence ID.
    fn mate_reference_sequence_id(&self, header: &Header) -> Option<&dyn ReferenceSequenceId>;

    /// Returns the mate alignment start.
    fn mate_alignment_start(&self) -> Option<&dyn Position>;

    /// Returns the template length.
    fn template_length(&self) -> &dyn TemplateLength;

    /// Returns the sequence.
    fn sequence(&self) -> &dyn Sequence;

    /// Returns the quality scores.
    fn quality_scores(&self) -> &dyn QualityScores;

    /// Returns the data.
    fn data(&self) -> &dyn Data;
}
