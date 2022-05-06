//! BAM record builder.

use noodles_core::Position;
use noodles_sam as sam;

use super::Record;

/// A BAM record builder.
#[derive(Debug)]
pub struct Builder {
    reference_sequence_id: Option<usize>,
    position: Option<Position>,
    mapping_quality: Option<sam::alignment::record::MappingQuality>,
    flags: sam::alignment::record::Flags,
    mate_reference_sequence_id: Option<usize>,
    mate_position: Option<Position>,
    template_length: i32,
    read_name: Option<sam::alignment::record::ReadName>,
    cigar: sam::alignment::record::Cigar,
    sequence: sam::alignment::record::Sequence,
    quality_scores: sam::alignment::record::QualityScores,
    data: sam::alignment::record::Data,
}

impl Builder {
    /// Sets a reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    ///
    /// let record = bam::Record::builder()
    ///     .set_reference_sequence_id(1)
    ///     .build();
    ///
    /// assert_eq!(record.reference_sequence_id(), Some(1));
    /// ```
    pub fn set_reference_sequence_id(mut self, reference_sequence_id: usize) -> Self {
        self.reference_sequence_id = Some(reference_sequence_id);
        self
    }

    /// Sets a position.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_core::Position;
    ///
    /// let record = bam::Record::builder()
    ///     .set_position(Position::try_from(8)?)
    ///     .build();
    ///
    /// assert_eq!(record.position(), Position::new(8));
    /// # Ok::<(), noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn set_position(mut self, position: Position) -> Self {
        self.position = Some(position);
        self
    }

    /// Sets a mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{alignment::record::MappingQuality, AlignmentRecord};
    ///
    /// let record = bam::Record::builder()
    ///     .set_mapping_quality(MappingQuality::try_from(34)?)
    ///     .build();
    ///
    /// assert_eq!(record.mapping_quality(), MappingQuality::new(34));
    /// # Ok::<_, noodles_sam::alignment::record::mapping_quality::ParseError>(())
    /// ```
    pub fn set_mapping_quality(
        mut self,
        mapping_quality: sam::alignment::record::MappingQuality,
    ) -> Self {
        self.mapping_quality = Some(mapping_quality);
        self
    }

    /// Sets SAM flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{alignment::record::Flags, AlignmentRecord};
    ///
    /// let record = bam::Record::builder()
    ///     .set_flags(Flags::PAIRED | Flags::READ_1)
    ///     .build();
    ///
    /// assert_eq!(record.flags(), Flags::PAIRED | Flags::READ_1);
    /// ```
    pub fn set_flags(mut self, flags: sam::alignment::record::Flags) -> Self {
        self.flags = flags;
        self
    }

    /// Sets a mate reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    ///
    /// let record = bam::Record::builder()
    ///     .set_mate_reference_sequence_id(1)
    ///     .build();
    ///
    /// assert_eq!(record.mate_reference_sequence_id(), Some(1));
    /// ```
    pub fn set_mate_reference_sequence_id(mut self, mate_reference_sequence_id: usize) -> Self {
        self.mate_reference_sequence_id = Some(mate_reference_sequence_id);
        self
    }

    /// Sets a mate position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_core::Position;
    ///
    /// let record = bam::Record::builder()
    ///     .set_mate_position(Position::try_from(13)?)
    ///     .build();
    ///
    /// assert_eq!(record.mate_position(), Position::new(13));
    /// # Ok::<(), noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn set_mate_position(mut self, mate_position: Position) -> Self {
        self.mate_position = Some(mate_position);
        self
    }

    /// Sets a template length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::AlignmentRecord;
    /// let record = bam::Record::builder().set_template_length(144).build();
    /// assert_eq!(record.template_length(), 144);
    /// ```
    pub fn set_template_length(mut self, template_length: i32) -> Self {
        self.template_length = template_length;
        self
    }

    /// Sets a read name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{alignment::record::ReadName, AlignmentRecord};
    ///
    /// let read_name: ReadName = "r0".parse()?;
    ///
    /// let record = bam::Record::builder()
    ///     .set_read_name(read_name.clone())
    ///     .build();
    ///
    /// assert_eq!(record.read_name(), Some(&read_name));
    /// # Ok::<(), noodles_sam::alignment::record::read_name::ParseError>(())
    /// ```
    pub fn set_read_name(mut self, read_name: sam::alignment::record::ReadName) -> Self {
        self.read_name = Some(read_name);
        self
    }

    /// Sets a CIGAR.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{alignment::record::Cigar, AlignmentRecord};
    ///
    /// let cigar: Cigar = "36M".parse()?;
    ///
    /// let record = bam::Record::builder()
    ///     .set_cigar(cigar.clone())
    ///     .build();
    ///
    /// assert_eq!(record.cigar(), &cigar);
    /// Ok::<_, noodles_sam::alignment::record::cigar::ParseError>(())
    /// ```
    pub fn set_cigar(mut self, cigar: sam::alignment::record::Cigar) -> Self {
        self.cigar = cigar;
        self
    }

    /// Sets a sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{alignment::record::Sequence, AlignmentRecord};
    ///
    /// let sequence: Sequence = "ACGT".parse()?;
    ///
    /// let record = bam::Record::builder()
    ///     .set_sequence(sequence.clone())
    ///     .build();
    ///
    /// assert_eq!(record.sequence(), &sequence);
    /// # Ok::<_, noodles_sam::alignment::record::sequence::ParseError>(())
    /// ```
    pub fn set_sequence(mut self, sequence: sam::alignment::record::Sequence) -> Self {
        self.sequence = sequence;
        self
    }

    /// Sets quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{alignment::record::QualityScores, AlignmentRecord};
    ///
    /// let quality_scores: QualityScores = "NDLS".parse()?;
    ///
    /// let record = bam::Record::builder()
    ///     .set_quality_scores(quality_scores.clone())
    ///     .build();
    ///
    /// assert_eq!(record.quality_scores(), &quality_scores);
    /// # Ok::<_, noodles_sam::alignment::record::quality_scores::ParseError>(())
    /// ```
    pub fn set_quality_scores(
        mut self,
        quality_scores: sam::alignment::record::QualityScores,
    ) -> Self {
        self.quality_scores = quality_scores;
        self
    }

    /// Sets data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{alignment::record::Data, AlignmentRecord};
    ///
    /// let data: Data = "NH:i:1".parse()?;
    ///
    /// let record = bam::Record::builder()
    ///     .set_data(data.clone())
    ///     .build();
    ///
    /// assert_eq!(record.data(), &data);
    /// # Ok::<_, noodles_sam::alignment::record::data::ParseError>(())
    /// ```
    pub fn set_data(mut self, data: sam::alignment::record::Data) -> Self {
        self.data = data;
        self
    }

    /// Builds a BAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::builder().build();
    /// ```
    pub fn build(self) -> Record {
        Record {
            reference_sequence_id: self.reference_sequence_id,
            position: self.position,
            mapping_quality: self.mapping_quality,
            flags: self.flags,
            mate_reference_sequence_id: self.mate_reference_sequence_id,
            mate_position: self.mate_position,
            template_length: self.template_length,
            read_name: self.read_name,
            cigar: self.cigar,
            sequence: self.sequence,
            quality_scores: self.quality_scores,
            data: self.data,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            reference_sequence_id: None,
            position: None,
            mapping_quality: None,
            flags: sam::alignment::record::Flags::UNMAPPED,
            mate_reference_sequence_id: None,
            mate_position: None,
            template_length: 0,
            read_name: None,
            cigar: sam::alignment::record::Cigar::default(),
            sequence: sam::alignment::record::Sequence::default(),
            quality_scores: sam::alignment::record::QualityScores::default(),
            data: sam::alignment::record::Data::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use sam::alignment::record::{AlignmentQualityScores, AlignmentSequence};

    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert!(builder.reference_sequence_id.is_none());
        assert!(builder.position.is_none());
        assert!(builder.mapping_quality.is_none());
        assert_eq!(builder.flags, sam::alignment::record::Flags::UNMAPPED);
        assert!(builder.mate_reference_sequence_id.is_none());
        assert!(builder.mate_position.is_none());
        assert_eq!(builder.template_length, 0);
        assert!(builder.read_name.is_none());
        assert!(builder.cigar.is_empty());
        assert!(builder.sequence.is_empty());
        assert!(builder.quality_scores.is_empty());
        assert!(builder.data.is_empty());
    }
}
