use super::{
    Cigar, Data, Flags, MappingQuality, Position, QualityScores, ReadName, Record,
    ReferenceSequenceName, Sequence,
};

/// A SAM record builder.
#[derive(Debug)]
pub struct Builder {
    read_name: ReadName,
    flags: Flags,
    reference_sequence_name: Option<ReferenceSequenceName>,
    position: Option<Position>,
    mapping_quality: MappingQuality,
    cigar: Cigar,
    mate_reference_sequence_name: Option<ReferenceSequenceName>,
    mate_position: Option<Position>,
    template_len: i32,
    sequence: Sequence,
    quality_scores: QualityScores,
    data: Data,
}

impl Builder {
    /// Creates a SAM record builder.
    ///
    /// Typically, [`sam::Record::builder`] is used instead of calling
    /// [`sam::record::Builder::new`].
    ///
    /// [`sam::Record::builder`]: struct.Record.html#method.builder
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let builder = sam::Record::builder();
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets a SAM record read name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let record = sam::Record::builder()
    ///     .set_read_name("r0".parse()?)
    ///     .build();
    ///
    /// assert_eq!(record.read_name().as_ref(), "r0");
    /// Ok::<(), sam::record::read_name::ParseError>(())
    /// ```
    pub fn set_read_name(mut self, read_name: ReadName) -> Self {
        self.read_name = read_name;
        self
    }

    /// Sets SAM record flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Flags};
    ///
    /// let record = sam::Record::builder()
    ///     .set_flags(Flags::PAIRED | Flags::READ_1)
    ///     .build();
    ///
    /// assert_eq!(record.flags(), Flags::PAIRED | Flags::READ_1);
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
    ///     .build();
    ///
    /// assert_eq!(record.reference_sequence_name().map(|name| name.as_str()), Some("sq0"));
    /// Ok::<(), sam::record::reference_sequence_name::ParseError>(())
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
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_sam::{self as sam, record::Position};
    ///
    /// let record = sam::Record::builder()
    ///     .set_position(Position::try_from(13)?)
    ///     .build();
    ///
    /// assert_eq!(record.position().map(i32::from), Some(13));
    /// # Ok::<(), sam::record::position::TryFromIntError>(())
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
    /// use noodles_sam::{self as sam, record::MappingQuality};
    ///
    /// let record = sam::Record::builder()
    ///     .set_mapping_quality(MappingQuality::from(34))
    ///     .build();
    ///
    /// assert_eq!(*record.mapping_quality(), Some(34));
    /// ```
    pub fn set_mapping_quality(mut self, mapping_quality: MappingQuality) -> Self {
        self.mapping_quality = mapping_quality;
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
    ///     .build();
    ///
    /// assert_eq!(**record.cigar(), [Op::new(Kind::Match, 36)]);
    /// Ok::<(), sam::record::cigar::ParseError>(())
    /// ```
    pub fn set_cigar(mut self, cigar: Cigar) -> Self {
        self.cigar = cigar;
        self
    }

    /// Sets a SAM record mate reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let record = sam::Record::builder()
    ///     .set_mate_reference_sequence_name("sq0".parse()?)
    ///     .build();
    ///
    /// assert_eq!(record.mate_reference_sequence_name().map(|name| name.as_str()), Some("sq0"));
    /// Ok::<(), sam::record::reference_sequence_name::ParseError>(())
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
    /// # use std::convert::TryFrom;
    /// use noodles_sam::{self as sam, record::Position};
    ///
    /// let record = sam::Record::builder()
    ///     .set_mate_position(Position::try_from(17)?)
    ///     .build();
    ///
    /// assert_eq!(record.mate_position().map(i32::from), Some(17));
    /// # Ok::<(), sam::record::position::TryFromIntError>(())
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
    /// use noodles_sam as sam;
    /// let record = sam::Record::builder().set_template_len(36).build();
    /// assert_eq!(record.template_len(), 36);
    /// ```
    pub fn set_template_len(mut self, template_len: i32) -> Self {
        self.template_len = template_len;
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
    ///     .set_sequence("ACGT".parse()?)
    ///     .build();
    ///
    /// assert_eq!(**record.sequence(), [Base::A, Base::C, Base::G, Base::T]);
    /// Ok::<(), sam::record::sequence::ParseError>(())
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
    /// # use std::convert::TryFrom;
    /// use noodles_sam::{self as sam, record::quality_scores::Score};
    ///
    /// let record = sam::Record::builder()
    ///     .set_quality_scores("NDLS".parse()?)
    ///     .build();
    ///
    /// assert_eq!(**record.quality_scores(), [
    ///     Score::try_from('N')?,
    ///     Score::try_from('D')?,
    ///     Score::try_from('L')?,
    ///     Score::try_from('S')?,
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
    /// use noodles_sam::{self as sam, record::data::{field::{Tag, Value}, Field}};
    ///
    /// let record = sam::Record::builder()
    ///     .set_data("NH:i:1\tRG:Z:rg0".parse()?)
    ///     .build();
    ///
    /// assert_eq!(**record.data(), [
    ///     Field::new(Tag::AlignmentHitCount, Value::Int32(1)),
    ///     Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
    /// ]);
    /// Ok::<(), sam::record::data::ParseError>(())
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
    /// let record = sam::Record::builder().build();
    /// ```
    pub fn build(self) -> Record {
        Record {
            qname: self.read_name,
            flag: self.flags,
            rname: self.reference_sequence_name,
            pos: self.position,
            mapq: self.mapping_quality,
            cigar: self.cigar,
            rnext: self.mate_reference_sequence_name,
            pnext: self.mate_position,
            tlen: self.template_len,
            seq: self.sequence,
            qual: self.quality_scores,
            data: self.data,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            read_name: ReadName::default(),
            flags: Flags::UNMAPPED,
            reference_sequence_name: Default::default(),
            position: Default::default(),
            mapping_quality: MappingQuality::default(),
            cigar: Cigar::default(),
            mate_reference_sequence_name: Default::default(),
            mate_position: Default::default(),
            template_len: Default::default(),
            sequence: Sequence::default(),
            quality_scores: QualityScores::default(),
            data: Data::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::convert::TryFrom;

    use crate::record::{cigar, data};

    use super::*;

    #[test]
    fn test_default() {
        let record = Builder::new().build();

        assert!(record.read_name().is_none());
        assert_eq!(record.flags(), Flags::UNMAPPED);
        assert!(record.reference_sequence_name().is_none());
        assert!(record.position().is_none());
        assert!(record.mapping_quality().is_none());
        assert!(record.cigar().is_empty());
        assert!(record.mate_reference_sequence_name().is_none());
        assert!(record.mate_position().is_none());
        assert_eq!(record.template_len(), 0);
        assert!(record.sequence().is_empty());
        assert!(record.quality_scores().is_empty());
        assert!(record.data().is_empty());
    }

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        let read_name: ReadName = "r0".parse()?;
        let reference_sequence_name: ReferenceSequenceName = "sq0".parse()?;
        let cigar = Cigar::from(vec![cigar::Op::new(cigar::op::Kind::Match, 4)]);
        let mate_reference_sequence_name = reference_sequence_name.clone();
        let sequence: Sequence = "ATCGATC".parse()?;
        let quality_scores: QualityScores = "NOODLES".parse()?;

        let data = Data::from(vec![data::Field::new(
            data::field::Tag::AlignmentHitCount,
            data::field::Value::Int32(1),
        )]);

        let record = Builder::new()
            .set_read_name(read_name.clone())
            .set_flags(Flags::PAIRED | Flags::READ_1)
            .set_reference_sequence_name(reference_sequence_name.clone())
            .set_position(Position::try_from(13)?)
            .set_mapping_quality(MappingQuality::from(37))
            .set_cigar(cigar)
            .set_mate_reference_sequence_name(mate_reference_sequence_name.clone())
            .set_mate_position(Position::try_from(17)?)
            .set_template_len(4)
            .set_sequence(sequence.clone())
            .set_quality_scores(quality_scores.clone())
            .set_data(data)
            .build();

        assert_eq!(record.read_name(), &read_name);
        assert_eq!(record.flags(), Flags::PAIRED | Flags::READ_1);
        assert_eq!(
            record.reference_sequence_name(),
            Some(&reference_sequence_name)
        );
        assert_eq!(record.position().map(i32::from), Some(13));
        assert_eq!(u8::from(record.mapping_quality()), 37);
        assert_eq!(record.cigar().len(), 1);

        assert_eq!(
            record.mate_reference_sequence_name(),
            Some(&mate_reference_sequence_name)
        );

        assert_eq!(record.mate_position().map(i32::from), Some(17));
        assert_eq!(record.template_len(), 4);
        assert_eq!(record.sequence(), &sequence);
        assert_eq!(record.quality_scores(), &quality_scores);
        assert_eq!(record.data().len(), 1);

        Ok(())
    }
}
