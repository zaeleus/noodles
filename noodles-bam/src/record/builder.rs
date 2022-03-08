//! BAM record builder.

use std::{error, fmt};

use noodles_sam as sam;

use super::{Cigar, Data, Record, ReferenceSequenceId};

/// A BAM record builder.
#[derive(Debug)]
pub struct Builder {
    ref_id: Option<ReferenceSequenceId>,
    pos: Option<sam::record::Position>,
    mapq: Option<sam::record::MappingQuality>,
    flag: sam::record::Flags,
    next_ref_id: Option<ReferenceSequenceId>,
    next_pos: Option<sam::record::Position>,
    tlen: i32,
    read_name: Option<sam::record::ReadName>,
    cigar: Cigar,
    seq: sam::record::Sequence,
    qual: sam::record::QualityScores,
    data: Data,
}

impl Builder {
    /// Sets a reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::{self as bam, record::ReferenceSequenceId};
    ///
    /// let record = bam::Record::builder()
    ///     .set_reference_sequence_id(ReferenceSequenceId::from(1))
    ///     .build()?;
    ///
    /// assert_eq!(
    ///     record.reference_sequence_id(),
    ///     Some(ReferenceSequenceId::from(1))
    /// );
    /// # Ok::<_, bam::record::builder::BuildError>(())
    /// ```
    pub fn set_reference_sequence_id(mut self, reference_sequence_id: ReferenceSequenceId) -> Self {
        self.ref_id = Some(reference_sequence_id);
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
    /// use noodles_sam::record::Position;
    ///
    /// let record = bam::Record::builder()
    ///     .set_position(Position::try_from(8)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.position().map(i32::from), Some(8));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_position(mut self, position: sam::record::Position) -> Self {
        self.pos = Some(position);
        self
    }

    /// Sets a mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{record::MappingQuality, AlignmentRecord};
    ///
    /// let record = bam::Record::builder()
    ///     .set_mapping_quality(MappingQuality::try_from(34)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.mapping_quality().map(u8::from), Some(34));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_mapping_quality(mut self, mapping_quality: sam::record::MappingQuality) -> Self {
        self.mapq = Some(mapping_quality);
        self
    }

    /// Sets SAM flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{record::Flags, AlignmentRecord};
    ///
    /// let record = bam::Record::builder()
    ///     .set_flags(Flags::PAIRED | Flags::READ_1)
    ///     .build()?;
    ///
    /// assert_eq!(record.flags(), Flags::PAIRED | Flags::READ_1);
    /// # Ok::<_, bam::record::builder::BuildError>(())
    /// ```
    pub fn set_flags(mut self, flags: sam::record::Flags) -> Self {
        self.flag = flags;
        self
    }

    /// Sets a mate reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::{self as bam, record::ReferenceSequenceId};
    ///
    /// let record = bam::Record::builder()
    ///     .set_mate_reference_sequence_id(ReferenceSequenceId::from(1))
    ///     .build()?;
    ///
    /// assert_eq!(
    ///     record.mate_reference_sequence_id(),
    ///     Some(ReferenceSequenceId::from(1))
    /// );
    /// # Ok::<_, bam::record::builder::BuildError>(())
    /// ```
    pub fn set_mate_reference_sequence_id(
        mut self,
        mate_reference_sequence_id: ReferenceSequenceId,
    ) -> Self {
        self.next_ref_id = Some(mate_reference_sequence_id);
        self
    }

    /// Sets a mate position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::record::Position;
    ///
    /// let record = bam::Record::builder()
    ///     .set_position(Position::try_from(13)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.position().map(i32::from), Some(13));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_mate_position(mut self, mate_position: sam::record::Position) -> Self {
        self.next_pos = Some(mate_position);
        self
    }

    /// Sets a template length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::AlignmentRecord;
    /// let record = bam::Record::builder().set_template_length(144).build()?;
    /// assert_eq!(record.template_length(), 144);
    /// # Ok::<_, bam::record::builder::BuildError>(())
    /// ```
    pub fn set_template_length(mut self, template_length: i32) -> Self {
        self.tlen = template_length;
        self
    }

    /// Sets a read name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{record::ReadName, AlignmentRecord};
    ///
    /// let read_name = ReadName::try_new("r0")?;
    ///
    /// let record = bam::Record::builder()
    ///     .set_read_name(read_name.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.read_name(), Some(&read_name));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_read_name(mut self, read_name: sam::record::ReadName) -> Self {
        self.read_name = Some(read_name);
        self
    }

    /// Sets a CIGAR.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::{self as bam, record::{cigar::Op, Cigar}};
    /// use noodles_sam::record::cigar::op::Kind;
    ///
    /// let cigar = Cigar::from(vec![Op::new(Kind::Match, 36)?]);
    ///
    /// let record = bam::Record::builder()
    ///     .set_cigar(cigar.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.cigar(), &cigar);
    /// Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_cigar(mut self, cigar: Cigar) -> Self {
        self.cigar = cigar;
        self
    }

    /// Sets a sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::record::Sequence;
    ///
    /// let sequence: Sequence = "ACGT".parse()?;
    ///
    /// let record = bam::Record::builder()
    ///     .set_sequence(sequence.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.sequence(), &sequence);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_sequence(mut self, sequence: sam::record::Sequence) -> Self {
        self.seq = sequence;
        self
    }

    /// Sets quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{record::{quality_scores::Score, QualityScores}, AlignmentRecord};
    ///
    /// let quality_scores: QualityScores = "NDLS".parse()?;
    ///
    /// let record = bam::Record::builder()
    ///     .set_quality_scores(quality_scores.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.quality_scores(), &quality_scores);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_quality_scores(mut self, quality_scores: sam::record::QualityScores) -> Self {
        self.qual = quality_scores;
        self
    }

    /// Sets data.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::{self as bam, record::Data};
    ///
    /// let data = Data::try_from(vec![
    ///     b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00, // NH:i:1
    /// ])?;
    ///
    /// let record = bam::Record::builder()
    ///     .set_data(data.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.data(), &data);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_data(mut self, data: Data) -> Self {
        self.data = data;
        self
    }

    /// Builds a BAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::builder().build()?;
    /// # Ok::<_, bam::record::builder::BuildError>(())
    /// ```
    pub fn build(self) -> Result<Record, BuildError> {
        Ok(Record {
            ref_id: self.ref_id,
            pos: self.pos,
            mapq: self.mapq,
            flag: self.flag,
            next_ref_id: self.next_ref_id,
            next_pos: self.next_pos,
            tlen: self.tlen,
            read_name: self.read_name,
            cigar: self.cigar,
            seq: self.seq,
            qual: self.qual,
            data: self.data,
        })
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            ref_id: None,
            pos: None,
            mapq: None,
            flag: sam::record::Flags::UNMAPPED,
            next_ref_id: None,
            next_pos: None,
            tlen: 0,
            read_name: None,
            cigar: Cigar::default(),
            seq: sam::record::Sequence::default(),
            qual: sam::record::QualityScores::default(),
            data: Data::default(),
        }
    }
}

/// An error returned when a BAM record fails to build.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BuildError {
    /// The CIGAR is invalid.
    InvalidCigar,
}

impl error::Error for BuildError {}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidCigar => f.write_str("invalid CIGAR"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert!(builder.ref_id.is_none());
        assert!(builder.pos.is_none());
        assert!(builder.mapq.is_none());
        assert_eq!(builder.flag, sam::record::Flags::UNMAPPED);
        assert!(builder.next_ref_id.is_none());
        assert!(builder.next_pos.is_none());
        assert_eq!(builder.tlen, 0);
        assert!(builder.read_name.is_none());
        assert!(builder.cigar.is_empty());
        assert!(builder.seq.is_empty());
        assert!(builder.qual.is_empty());
        assert!(builder.data.is_empty());
    }
}
