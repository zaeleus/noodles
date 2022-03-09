//! SAM record builder.

use std::{error, fmt};

use super::{
    Cigar, Data, Flags, MappingQuality, Position, QualityScores, ReadName, Record,
    ReferenceSequenceName, Sequence,
};

/// An error returned when a SAM record fails to build.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BuildError {
    /// The sequence length does not match the CIGAR string read length.
    SequenceLengthMismatch(usize, usize),
    /// The quality scores length does not match the sequence length.
    QualityScoresLengthMismatch(usize, usize),
}

impl error::Error for BuildError {}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SequenceLengthMismatch(sequence_len, cigar_read_len) => write!(
                f,
                "sequence length mismatch: expected {}, got {}",
                cigar_read_len, sequence_len
            ),
            Self::QualityScoresLengthMismatch(quality_scores_len, sequence_len) => write!(
                f,
                "quality scores length mismatch: expected {}, got {}",
                sequence_len, quality_scores_len
            ),
        }
    }
}

/// A SAM record builder.
#[derive(Debug)]
pub struct Builder {
    read_name: Option<ReadName>,
    flags: Flags,
    reference_sequence_name: Option<ReferenceSequenceName>,
    position: Option<Position>,
    mapping_quality: Option<MappingQuality>,
    cigar: Cigar,
    mate_reference_sequence_name: Option<ReferenceSequenceName>,
    mate_position: Option<Position>,
    template_length: i32,
    sequence: Sequence,
    quality_scores: QualityScores,
    data: Data,
}

impl Builder {
    /// Creates a SAM record builder.
    ///
    /// Typically, [`Record::builder`] is used instead of calling this.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let builder = sam::Record::builder();
    /// ```
    #[deprecated(since = "0.9.0", note = "Use `Record::builder` instead.")]
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets a SAM record read name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, AlignmentRecord};
    ///
    /// let record = sam::Record::builder()
    ///     .set_read_name("r0".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(record.read_name().map(|name| name.as_ref()), Some("r0"));
    /// Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_read_name(mut self, read_name: ReadName) -> Self {
        self.read_name = Some(read_name);
        self
    }

    /// Sets SAM record flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Flags, AlignmentRecord};
    ///
    /// let record = sam::Record::builder()
    ///     .set_flags(Flags::PAIRED | Flags::READ_1)
    ///     .build()?;
    ///
    /// assert_eq!(record.flags(), Flags::PAIRED | Flags::READ_1);
    /// # Ok::<(), sam::record::builder::BuildError>(())
    /// ```
    pub fn set_flags(mut self, flags: Flags) -> Self {
        self.flags = flags;
        self
    }

    /// Sets a SAM record reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let record = sam::Record::builder()
    ///     .set_reference_sequence_name("sq0".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(record.reference_sequence_name().map(|name| name.as_str()), Some("sq0"));
    /// Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_reference_sequence_name(
        mut self,
        reference_sequence_name: ReferenceSequenceName,
    ) -> Self {
        self.reference_sequence_name = Some(reference_sequence_name);
        self
    }

    /// Sets a SAM record position.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Position};
    ///
    /// let record = sam::Record::builder()
    ///     .set_position(Position::try_from(13)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.position().map(i32::from), Some(13));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_position(mut self, position: Position) -> Self {
        self.position = Some(position);
        self
    }

    /// Sets a SAM record mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::MappingQuality, AlignmentRecord};
    ///
    /// let record = sam::Record::builder()
    ///     .set_mapping_quality(MappingQuality::try_from(34)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.mapping_quality().map(u8::from), Some(34));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_mapping_quality(mut self, mapping_quality: MappingQuality) -> Self {
        self.mapping_quality = Some(mapping_quality);
        self
    }

    /// Sets a SAM record CIGAR.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::cigar::{op::Kind, Op}};
    ///
    /// let record = sam::Record::builder()
    ///     .set_cigar("36M".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(**record.cigar(), [Op::new(Kind::Match, 36)]);
    /// Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_cigar(mut self, cigar: Cigar) -> Self {
        self.cigar = cigar;
        self
    }

    /// Sets a SAM record mate reference sequence name.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let record = sam::Record::builder()
    ///     .set_mate_reference_sequence_name("sq0".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(record.mate_reference_sequence_name().map(|name| name.as_str()), Some("sq0"));
    /// Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_mate_reference_sequence_name(
        mut self,
        mate_reference_sequence_name: ReferenceSequenceName,
    ) -> Self {
        self.mate_reference_sequence_name = Some(mate_reference_sequence_name);
        self
    }

    /// Sets a SAM record mate position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Position};
    ///
    /// let record = sam::Record::builder()
    ///     .set_mate_position(Position::try_from(17)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.mate_position().map(i32::from), Some(17));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_mate_position(mut self, mate_position: Position) -> Self {
        self.mate_position = Some(mate_position);
        self
    }

    /// Sets a SAM record template length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, AlignmentRecord};
    /// let record = sam::Record::builder().set_template_length(36).build()?;
    /// assert_eq!(record.template_length(), 36);
    /// # Ok::<(), sam::record::builder::BuildError>(())
    /// ```
    pub fn set_template_length(mut self, template_length: i32) -> Self {
        self.template_length = template_length;
        self
    }

    /// Sets a SAM record sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::sequence::Base};
    ///
    /// let record = sam::Record::builder()
    ///     .set_cigar("4M".parse()?)
    ///     .set_sequence("ACGT".parse()?)
    ///     .set_quality_scores("NDLS".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(**record.sequence(), [Base::A, Base::C, Base::G, Base::T]);
    /// Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_sequence(mut self, sequence: Sequence) -> Self {
        self.sequence = sequence;
        self
    }

    /// Sets SAM record quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::quality_scores::Score, AlignmentRecord};
    ///
    /// let record = sam::Record::builder()
    ///     .set_cigar("2M".parse()?)
    ///     .set_sequence("AC".parse()?)
    ///     .set_quality_scores("ND".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(**record.quality_scores(), [
    ///     Score::try_from('N')?,
    ///     Score::try_from('D')?,
    /// ]);
    /// Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_quality_scores(mut self, quality_scores: QualityScores) -> Self {
        self.quality_scores = quality_scores;
        self
    }

    /// Sets SAM record data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::{data, Data}};
    ///
    /// let record = sam::Record::default();
    /// assert!(record.data().is_empty());
    ///
    /// let data = Data::try_from(vec![data::Field::new(
    ///     data::field::Tag::AlignmentHitCount,
    ///     data::field::Value::Int(1),
    /// )])?;
    /// let record = sam::Record::builder().set_data(data).build()?;
    /// assert_eq!(record.data().to_string(), "NH:i:1");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_data(mut self, data: Data) -> Self {
        self.data = data;
        self
    }

    /// Builds a SAM record.
    ///
    /// # Example
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::Record::builder().build()?;
    /// # Ok::<(), sam::record::builder::BuildError>(())
    /// ```
    pub fn build(self) -> Result<Record, BuildError> {
        // § 1.4 The alignment section: mandatory fields (2021-06-03): "If not a '*', the length of
        // the sequence must equal the sum of lengths of `M/I/S/=/X` operations in `CIGAR`."
        if !self.flags.is_unmapped() && !self.sequence.is_empty() {
            let sequence_len = self.sequence.len();
            let cigar_read_len = self.cigar.read_len();

            if sequence_len != cigar_read_len {
                return Err(BuildError::SequenceLengthMismatch(
                    sequence_len,
                    cigar_read_len,
                ));
            }
        }

        // § 1.4 The alignment section: mandatory fields (2021-06-03): "If not a '*', `SEQ` must
        // not be a '*' and the length of the quality string ought to equal the length of `SEQ`."
        if !self.quality_scores.is_empty() {
            let quality_scores_len = self.quality_scores.len();
            let sequence_len = self.sequence.len();

            if quality_scores_len != sequence_len {
                return Err(BuildError::QualityScoresLengthMismatch(
                    quality_scores_len,
                    sequence_len,
                ));
            }
        }

        Ok(Record {
            read_name: self.read_name,
            flags: self.flags,
            reference_sequence_name: self.reference_sequence_name,
            position: self.position,
            mapping_quality: self.mapping_quality,
            cigar: self.cigar,
            mate_reference_sequence_name: self.mate_reference_sequence_name,
            mate_position: self.mate_position,
            template_length: self.template_length,
            sequence: self.sequence,
            quality_scores: self.quality_scores,
            data: self.data,
        })
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            read_name: None,
            flags: Flags::UNMAPPED,
            reference_sequence_name: None,
            position: None,
            mapping_quality: None,
            cigar: Cigar::default(),
            mate_reference_sequence_name: None,
            mate_position: None,
            template_length: 0,
            sequence: Sequence::default(),
            quality_scores: QualityScores::default(),
            data: Data::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        record::{cigar, data},
        AlignmentRecord,
    };

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert!(builder.read_name.is_none());
        assert_eq!(builder.flags, Flags::UNMAPPED);
        assert!(builder.reference_sequence_name.is_none());
        assert!(builder.position.is_none());
        assert!(builder.mapping_quality.is_none());
        assert!(builder.cigar.is_empty());
        assert!(builder.mate_reference_sequence_name.is_none());
        assert!(builder.mate_position.is_none());
        assert_eq!(builder.template_length, 0);
        assert!(builder.sequence.is_empty());
        assert!(builder.quality_scores.is_empty());
        assert!(builder.data.is_empty());
    }

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        let read_name: ReadName = "r0".parse()?;
        let reference_sequence_name: ReferenceSequenceName = "sq0".parse()?;
        let cigar = Cigar::from(vec![cigar::Op::new(cigar::op::Kind::Match, 4)]);
        let mate_reference_sequence_name = reference_sequence_name.clone();
        let sequence: Sequence = "ATCG".parse()?;
        let quality_scores: QualityScores = "NDLS".parse()?;

        let data = Data::try_from(vec![data::Field::new(
            data::field::Tag::AlignmentHitCount,
            data::field::Value::Int(1),
        )])?;

        let record = Builder::default()
            .set_read_name(read_name.clone())
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .set_reference_sequence_name(reference_sequence_name.clone())
            .set_position(Position::try_from(13)?)
            .set_mapping_quality(MappingQuality::try_from(37)?)
            .set_cigar(cigar)
            .set_mate_reference_sequence_name(mate_reference_sequence_name.clone())
            .set_mate_position(Position::try_from(17)?)
            .set_template_length(4)
            .set_sequence(sequence.clone())
            .set_quality_scores(quality_scores.clone())
            .set_data(data)
            .build()?;

        assert_eq!(record.read_name(), Some(&read_name));
        assert_eq!(record.flags(), Flags::SEGMENTED | Flags::FIRST_SEGMENT);
        assert_eq!(
            record.reference_sequence_name(),
            Some(&reference_sequence_name)
        );
        assert_eq!(record.position().map(i32::from), Some(13));
        assert_eq!(record.mapping_quality().map(u8::from), Some(37));
        assert_eq!(record.cigar().len(), 1);

        assert_eq!(
            record.mate_reference_sequence_name(),
            Some(&mate_reference_sequence_name)
        );

        assert_eq!(record.mate_position().map(i32::from), Some(17));
        assert_eq!(record.template_length(), 4);
        assert_eq!(record.sequence(), &sequence);
        assert_eq!(record.quality_scores(), &quality_scores);
        assert_eq!(record.data().len(), 1);

        Ok(())
    }

    #[test]
    fn test_build_with_quality_scores_length_mismatch() -> Result<(), Box<dyn std::error::Error>> {
        let builder = Builder::default()
            .set_cigar("4M".parse()?)
            .set_sequence("ACGT".parse()?)
            .set_quality_scores("NDL".parse()?);

        assert_eq!(
            builder.build(),
            Err(BuildError::QualityScoresLengthMismatch(3, 4))
        );

        Ok(())
    }
}
