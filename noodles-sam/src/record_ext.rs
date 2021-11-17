use std::io;

use super::{
    header::{ReferenceSequence, ReferenceSequences},
    record::Position,
};

/// SAM(-like) record extensions.
pub trait RecordExt {
    /// Returns the associated reference sequence.
    fn reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>>;

    /// Returns the start position.
    fn alignment_start(&self) -> Option<Position>;

    /// Calculates the alignment span over the reference sequence.
    fn alignment_span(&self) -> io::Result<u32>;

    /// Calculates the end position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     self as sam,
    ///     record::{cigar::{op::Kind, Op}, Cigar, Position},
    ///     RecordExt,
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

        let len = match self.alignment_span().and_then(|len| {
            i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        }) {
            Ok(len) => len,
            Err(e) => return Some(Err(e)),
        };

        let end = start + len - 1;

        Some(Position::try_from(end).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
    }

    /// Returns the associated reference sequence of the mate.
    fn mate_reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>>;
}
