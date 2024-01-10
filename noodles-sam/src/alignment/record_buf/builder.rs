use noodles_core::Position;

use super::{Cigar, Data, Name, QualityScores, RecordBuf, Sequence};
use crate::record::{Flags, MappingQuality};

/// An alignment record builder.
#[derive(Debug)]
pub struct Builder {
    name: Option<Name>,
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

impl Builder {
    /// Sets the read name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, alignment::record_buf::Name};
    ///
    /// let name = Name::from(b"r1");
    ///
    /// let record = sam::alignment::RecordBuf::builder()
    ///     .set_name(name.clone())
    ///     .build();
    ///
    /// assert_eq!(record.name(), Some(&name));
    /// ```
    pub fn set_name(mut self, name: Name) -> Self {
        self.name = Some(name);
        self
    }

    /// Sets the flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Flags};
    ///
    /// let record = sam::alignment::RecordBuf::builder()
    ///     .set_flags(Flags::empty())
    ///     .build();
    ///
    /// assert_eq!(record.flags(), Flags::empty());
    /// ```
    pub fn set_flags(mut self, flags: Flags) -> Self {
        self.flags = flags;
        self
    }

    /// Sets the reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let record = sam::alignment::RecordBuf::builder()
    ///     .set_reference_sequence_id(0)
    ///     .build();
    ///
    /// assert_eq!(record.reference_sequence_id(), Some(0));
    /// ```
    pub fn set_reference_sequence_id(mut self, reference_sequence_id: usize) -> Self {
        self.reference_sequence_id = Some(reference_sequence_id);
        self
    }

    /// Sets the alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam as sam;
    ///
    /// let record = sam::alignment::RecordBuf::builder()
    ///     .set_alignment_start(Position::MIN)
    ///     .build();
    ///
    /// assert_eq!(record.alignment_start(), Some(Position::MIN));
    /// ```
    pub fn set_alignment_start(mut self, alignment_start: Position) -> Self {
        self.alignment_start = Some(alignment_start);
        self
    }

    /// Sets the mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::MappingQuality};
    ///
    /// let record = sam::alignment::RecordBuf::builder()
    ///     .set_mapping_quality(MappingQuality::MIN)
    ///     .build();
    ///
    /// assert_eq!(record.mapping_quality(), Some(MappingQuality::MIN));
    /// ```
    pub fn set_mapping_quality(mut self, mapping_quality: MappingQuality) -> Self {
        self.mapping_quality = Some(mapping_quality);
        self
    }

    /// Sets the CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     self as sam,
    ///     alignment::{
    ///         record::cigar::{op::Kind, Op},
    ///         record_buf::Cigar,
    ///     },
    /// };
    ///
    /// let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
    ///
    /// let record = sam::alignment::RecordBuf::builder()
    ///     .set_cigar(cigar.clone())
    ///     .build();
    ///
    /// assert_eq!(record.cigar(), &cigar);
    /// ```
    pub fn set_cigar(mut self, cigar: Cigar) -> Self {
        self.cigar = cigar;
        self
    }

    /// Sets the mate reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let record = sam::alignment::RecordBuf::builder()
    ///     .set_mate_reference_sequence_id(0)
    ///     .build();
    ///
    /// assert_eq!(record.mate_reference_sequence_id(), Some(0));
    /// ```
    pub fn set_mate_reference_sequence_id(mut self, mate_reference_sequence_id: usize) -> Self {
        self.mate_reference_sequence_id = Some(mate_reference_sequence_id);
        self
    }

    /// Sets the mate alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam as sam;
    ///
    /// let record = sam::alignment::RecordBuf::builder()
    ///     .set_mate_alignment_start(Position::MIN)
    ///     .build();
    ///
    /// assert_eq!(record.mate_alignment_start(), Some(Position::MIN));
    /// ```
    pub fn set_mate_alignment_start(mut self, mate_alignment_start: Position) -> Self {
        self.mate_alignment_start = Some(mate_alignment_start);
        self
    }

    /// Sets the template length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let record = sam::alignment::RecordBuf::builder()
    ///     .set_template_length(4)
    ///     .build();
    ///
    /// assert_eq!(record.template_length(), 4);
    /// ```
    pub fn set_template_length(mut self, template_length: i32) -> Self {
        self.template_length = template_length;
        self
    }

    /// Sets the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, alignment::record_buf::Sequence};
    ///
    /// let sequence = Sequence::from(b"ACGT".to_vec());
    ///
    /// let record = sam::alignment::RecordBuf::builder()
    ///     .set_sequence(sequence.clone())
    ///     .build();
    ///
    /// assert_eq!(record.sequence(), &sequence);
    /// ```
    pub fn set_sequence(mut self, sequence: Sequence) -> Self {
        self.sequence = sequence;
        self
    }

    /// Sets the quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, alignment::record_buf::QualityScores};
    ///
    /// let quality_scores = QualityScores::from(vec![45, 35, 43, 50]);
    ///
    /// let record = sam::alignment::RecordBuf::builder()
    ///     .set_quality_scores(quality_scores.clone())
    ///     .build();
    ///
    /// assert_eq!(record.quality_scores(), &quality_scores);
    /// ```
    pub fn set_quality_scores(mut self, quality_scores: QualityScores) -> Self {
        self.quality_scores = quality_scores;
        self
    }

    /// Sets the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     self as sam,
    ///     alignment::{
    ///         record::data::field::tag,
    ///         record_buf::{data::field::Value, Data},
    ///     },
    /// };
    ///
    /// let data: Data = [(tag::ALIGNMENT_HIT_COUNT, Value::from(1))]
    ///     .into_iter()
    ///     .collect();
    ///
    /// let record = sam::alignment::RecordBuf::builder()
    ///     .set_data(data.clone())
    ///     .build();
    ///
    /// assert_eq!(record.data(), &data);
    /// ```
    pub fn set_data(mut self, data: Data) -> Self {
        self.data = data;
        self
    }

    /// Builds the alignment record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::RecordBuf::builder().build();
    /// assert_eq!(record, sam::alignment::RecordBuf::default());
    /// ```
    pub fn build(self) -> RecordBuf {
        RecordBuf {
            name: self.name,
            flags: self.flags,
            reference_sequence_id: self.reference_sequence_id,
            alignment_start: self.alignment_start,
            mapping_quality: self.mapping_quality,
            cigar: self.cigar,
            mate_reference_sequence_id: self.mate_reference_sequence_id,
            mate_alignment_start: self.mate_alignment_start,
            template_length: self.template_length,
            sequence: self.sequence,
            quality_scores: self.quality_scores,
            data: self.data,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            name: None,
            flags: Flags::UNMAPPED,
            reference_sequence_id: None,
            alignment_start: None,
            mapping_quality: None,
            cigar: Cigar::default(),
            mate_reference_sequence_id: None,
            mate_alignment_start: None,
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

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert!(builder.name.is_none());
        assert_eq!(builder.flags, Flags::UNMAPPED);
        assert!(builder.reference_sequence_id.is_none());
        assert!(builder.alignment_start.is_none());
        assert!(builder.mapping_quality.is_none());
        assert!(builder.cigar.as_ref().is_empty());
        assert!(builder.mate_reference_sequence_id.is_none());
        assert!(builder.mate_alignment_start.is_none());
        assert_eq!(builder.template_length, 0);
        assert!(builder.sequence.is_empty());
        assert!(builder.quality_scores.is_empty());
        assert!(builder.data.is_empty());
    }
}
