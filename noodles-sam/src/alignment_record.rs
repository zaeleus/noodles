use std::io;

use noodles_core::Position;

use super::{
    header::{ReferenceSequence, ReferenceSequences},
    record::{Cigar, Flags, MappingQuality, QualityScores, ReadName, Sequence},
};

/// An alignment record.
pub trait AlignmentRecord {
    /// Returns the read name.
    ///
    /// This is also called the query name.
    fn read_name(&self) -> Option<&ReadName>;

    /// Returns the associated reference sequence.
    fn reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>>;

    /// Returns the flags.
    fn flags(&self) -> Flags;

    /// Returns the start position.
    fn alignment_start(&self) -> Option<Position>;

    /// Calculates the alignment span over the reference sequence.
    fn alignment_span(&self) -> usize;

    /// Calculates the end position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam::{self as sam, AlignmentRecord};
    ///
    /// let record = sam::Record::builder()
    ///     .set_position(Position::try_from(8)?)
    ///     .set_cigar("5M".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(record.alignment_end(), Position::new(12));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    fn alignment_end(&self) -> Option<Position> {
        self.alignment_start().and_then(|alignment_start| {
            let end = usize::from(alignment_start) + self.alignment_span() - 1;
            Position::new(end)
        })
    }

    /// Returns the mapping quality.
    fn mapping_quality(&self) -> Option<MappingQuality>;

    /// Returns the CIGAR operations.
    fn cigar(&self) -> &Cigar;

    /// Returns the associated reference sequence of the mate.
    fn mate_reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>>;

    /// Returns the start position of the mate.
    fn mate_alignment_start(&self) -> Option<Position>;

    /// Returns the template length.
    fn template_length(&self) -> i32;

    /// Returns the sequence.
    fn sequence(&self) -> &Sequence;

    /// Returns the quality scores.
    fn quality_scores(&self) -> &QualityScores;
}
