use super::record::{
    Cigar, Data, Flags, MappingQuality, Position, QualityScores, ReadName, ReferenceSequenceId,
    Sequence, TemplateLength,
};
use crate::Header;

/// An alignment record.
pub trait AnyRecord {
    /// Returns the read name.
    fn read_name(&self) -> Option<Box<dyn ReadName + '_>>;

    /// Returns the flags.
    fn flags(&self) -> Box<dyn Flags + '_>;

    /// Returns the reference sequence ID.
    fn reference_sequence_id(&self, header: &Header) -> Option<Box<dyn ReferenceSequenceId + '_>>;

    /// Returns the alignment start.
    fn alignment_start(&self) -> Option<Box<dyn Position + '_>>;

    /// Returns the mapping quality.
    fn mapping_quality(&self) -> Option<Box<dyn MappingQuality + '_>>;

    /// Returns the CIGAR operations.
    fn cigar(&self, header: &Header) -> Box<dyn Cigar + '_>;

    /// Returns the mate reference sequence ID.
    fn mate_reference_sequence_id(
        &self,
        header: &Header,
    ) -> Option<Box<dyn ReferenceSequenceId + '_>>;

    /// Returns the mate alignment start.
    fn mate_alignment_start(&self) -> Option<Box<dyn Position + '_>>;

    /// Returns the template length.
    fn template_length(&self) -> Box<dyn TemplateLength + '_>;

    /// Returns the sequence.
    fn sequence(&self) -> Box<dyn Sequence + '_>;

    /// Returns the quality scores.
    fn quality_scores(&self) -> Box<dyn QualityScores + '_>;

    /// Returns the data.
    fn data(&self) -> Box<dyn Data + '_>;
}
