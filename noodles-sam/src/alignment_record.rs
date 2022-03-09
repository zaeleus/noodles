use std::io;

use super::{
    header::{ReferenceSequence, ReferenceSequences},
    record::{Flags, MappingQuality, Position, QualityScores, ReadName},
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
    /// use noodles_sam::{
    ///     self as sam,
    ///     record::{cigar::{op::Kind, Op}, Cigar, Position},
    ///     AlignmentRecord,
    /// };
    ///
    /// let record = sam::Record::builder()
    ///     .set_position(Position::try_from(8)?)
    ///     .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 5)]))
    ///     .build()?;
    ///
    /// let actual = record.alignment_end().transpose()?;
    /// let expected = Position::try_from(12).map(Some)?;
    /// assert_eq!(actual, expected);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    fn alignment_end(&self) -> Option<io::Result<Position>> {
        let start = self.alignment_start().map(i32::from)?;
        let len = self.alignment_span() as i32;
        let end = start + len - 1;
        Some(Position::try_from(end).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
    }

    /// Returns the mapping quality.
    fn mapping_quality(&self) -> Option<MappingQuality>;

    /// Returns the associated reference sequence of the mate.
    fn mate_reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>>;

    /// Returns the start position of the mate.
    fn mate_alignment_start(&self) -> Option<Position>;

    /// Returns the template length.
    fn template_length(&self) -> i32;

    /// Returns the quality scores.
    fn quality_scores(&self) -> &QualityScores;
}
