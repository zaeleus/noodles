use std::io;

use noodles_core::Position;

use super::{
    alignment::record::{
        AlignmentQualityScores, AlignmentSequence, Cigar, Data, Flags, MappingQuality, ReadName,
    },
    header::{ReferenceSequence, ReferenceSequences},
};

/// An alignment record.
pub trait AlignmentRecord {
    /// The sequence returned.
    type Sequence: AlignmentSequence;

    /// The quality scores returned.
    type QualityScores: AlignmentQualityScores;

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
    ///     .build();
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
    fn sequence(&self) -> &Self::Sequence;

    /// Returns the quality scores.
    fn quality_scores(&self) -> &Self::QualityScores;

    /// Returns the data fields.
    fn data(&self) -> &Data;
}

impl<R> AlignmentRecord for Box<R>
where
    R: AlignmentRecord + ?Sized,
{
    type Sequence = R::Sequence;
    type QualityScores = R::QualityScores;

    fn read_name(&self) -> Option<&ReadName> {
        (**self).read_name()
    }

    fn reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>> {
        (**self).reference_sequence(reference_sequences)
    }

    fn flags(&self) -> Flags {
        (**self).flags()
    }

    fn alignment_start(&self) -> Option<Position> {
        (**self).alignment_start()
    }

    fn alignment_span(&self) -> usize {
        (**self).alignment_span()
    }

    fn mapping_quality(&self) -> Option<MappingQuality> {
        (**self).mapping_quality()
    }

    fn cigar(&self) -> &Cigar {
        (**self).cigar()
    }

    fn mate_reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>> {
        (**self).mate_reference_sequence(reference_sequences)
    }

    fn mate_alignment_start(&self) -> Option<Position> {
        (**self).mate_alignment_start()
    }

    fn template_length(&self) -> i32 {
        (**self).template_length()
    }

    fn sequence(&self) -> &Self::Sequence {
        (**self).sequence()
    }

    fn quality_scores(&self) -> &Self::QualityScores {
        (**self).quality_scores()
    }

    fn data(&self) -> &Data {
        (**self).data()
    }
}

/// Any alignment record.
pub trait AnyAlignmentRecord {
    /// Returns the read name.
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
    fn alignment_end(&self) -> Option<Position>;

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
    fn sequence(&self) -> &dyn AlignmentSequence;

    /// Returns the quality scores.
    fn quality_scores(&self) -> &dyn AlignmentQualityScores;

    /// Returns the data fields.
    fn data(&self) -> &Data;
}

impl<T> AnyAlignmentRecord for T
where
    T: AlignmentRecord,
{
    fn read_name(&self) -> Option<&ReadName> {
        self.read_name()
    }

    fn reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>> {
        self.reference_sequence(reference_sequences)
    }

    fn flags(&self) -> Flags {
        self.flags()
    }

    fn alignment_start(&self) -> Option<Position> {
        self.alignment_start()
    }

    fn alignment_span(&self) -> usize {
        self.alignment_span()
    }

    fn alignment_end(&self) -> Option<Position> {
        self.alignment_end()
    }

    fn mapping_quality(&self) -> Option<MappingQuality> {
        self.mapping_quality()
    }

    fn cigar(&self) -> &Cigar {
        self.cigar()
    }

    fn mate_reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>> {
        self.mate_reference_sequence(reference_sequences)
    }

    fn mate_alignment_start(&self) -> Option<Position> {
        self.mate_alignment_start()
    }

    fn template_length(&self) -> i32 {
        self.template_length()
    }

    fn sequence(&self) -> &dyn AlignmentSequence {
        self.sequence()
    }

    fn quality_scores(&self) -> &dyn AlignmentQualityScores {
        self.quality_scores()
    }

    fn data(&self) -> &Data {
        self.data()
    }
}
