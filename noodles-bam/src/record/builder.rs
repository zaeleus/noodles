use noodles_sam as sam;

use super::{Cigar, Data, QualityScores, Record, ReferenceSequenceId, Sequence};

/// A BAM record builder.
#[derive(Debug)]
pub struct Builder {
    ref_id: Option<ReferenceSequenceId>,
    pos: Option<sam::record::Position>,
    mapq: sam::record::MappingQuality,
    flag: sam::record::Flags,
    next_ref_id: Option<ReferenceSequenceId>,
    next_pos: Option<sam::record::Position>,
    tlen: i32,
    read_name: Vec<u8>,
    cigar: Cigar,
    seq: Sequence,
    qual: QualityScores,
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
    ///     .set_reference_sequence_id(ReferenceSequenceId::try_from(1)?)
    ///     .build();
    ///
    /// assert_eq!(record.reference_sequence_id().map(i32::from), Some(1));
    /// # Ok::<_, bam::record::reference_sequence_id::TryFromIntError>(())
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
    ///     .build();
    ///
    /// assert_eq!(record.position().map(i32::from), Some(8));
    /// # Ok::<_, noodles_sam::record::position::TryFromIntError>(())
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
    /// use noodles_sam::record::MappingQuality;
    ///
    /// let record = bam::Record::builder()
    ///     .set_mapping_quality(MappingQuality::from(34))
    ///     .build();
    ///
    /// assert_eq!(*record.mapping_quality(), Some(34));
    /// ```
    pub fn set_mapping_quality(mut self, mapping_quality: sam::record::MappingQuality) -> Self {
        self.mapq = mapping_quality;
        self
    }

    /// Sets SAM flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::record::Flags;
    ///
    /// let record = bam::Record::builder()
    ///     .set_flags(Flags::PAIRED | Flags::READ_1)
    ///     .build();
    ///
    /// assert_eq!(record.flags(), Flags::PAIRED | Flags::READ_1);
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
    ///     .set_mate_reference_sequence_id(ReferenceSequenceId::try_from(1)?)
    ///     .build();
    ///
    /// assert_eq!(record.mate_reference_sequence_id().map(i32::from), Some(1));
    /// # Ok::<_, bam::record::reference_sequence_id::TryFromIntError>(())
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
    ///     .build();
    ///
    /// assert_eq!(record.position().map(i32::from), Some(13));
    /// # Ok::<_, noodles_sam::record::position::TryFromIntError>(())
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
    /// let record = bam::Record::builder().set_template_length(144).build();
    /// assert_eq!(record.template_length(), 144);
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
    /// # use std::ffi;
    /// use noodles_bam as bam;
    ///
    /// let record = bam::Record::builder()
    ///     .set_read_name(b"r0\x00".to_vec())
    ///     .build();
    ///
    /// assert_eq!(record.read_name()?.to_bytes(), b"r0");
    /// # Ok::<(), ffi::FromBytesWithNulError>(())
    /// ```
    pub fn set_read_name(mut self, read_name: Vec<u8>) -> Self {
        self.read_name = read_name;
        self
    }

    /// Sets a CIGAR.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::{self as bam, record::Cigar};
    ///
    /// let cigar = Cigar::from(vec![0x00000240]); // 36M
    ///
    /// let record = bam::Record::builder()
    ///     .set_cigar(cigar.clone())
    ///     .build();
    ///
    /// assert_eq!(record.cigar(), &cigar);
    /// Ok::<_, bam::record::cigar::op::LengthError>(())
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
    /// use noodles_bam::{self as bam, record::Sequence};
    ///
    /// let sequence = Sequence::new(vec![0x12, 0x48], 4); // A
    ///
    /// let record = bam::Record::builder()
    ///     .set_sequence(sequence.clone())
    ///     .build();
    ///
    /// assert_eq!(record.sequence(), &sequence);
    /// ```
    pub fn set_sequence(mut self, sequence: Sequence) -> Self {
        self.seq = sequence;
        self
    }

    /// Sets quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::{self as bam, record::QualityScores};
    ///
    /// let quality_scores = QualityScores::from(vec![45, 35, 43, 50]); // NDLS
    ///
    /// let record = bam::Record::builder()
    ///     .set_quality_scores(quality_scores.clone())
    ///     .build();
    ///
    /// assert_eq!(record.quality_scores(), &quality_scores);
    /// ```
    pub fn set_quality_scores(mut self, quality_scores: QualityScores) -> Self {
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
    ///     .build();
    ///
    /// assert_eq!(record.data(), &data);
    /// # Ok::<_, io::Error>(())
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
    /// let record = bam::Record::builder().build();
    /// ```
    pub fn build(self) -> Record {
        use super::{reference_sequence_id, UNMAPPED_POSITION};

        let ref_id = self
            .ref_id
            .map(i32::from)
            .unwrap_or(reference_sequence_id::UNMAPPED);

        let pos = self
            .pos
            .map(|p| i32::from(p) - 1)
            .unwrap_or(UNMAPPED_POSITION);

        let next_ref_id = self
            .next_ref_id
            .map(i32::from)
            .unwrap_or(reference_sequence_id::UNMAPPED);

        let next_pos = self
            .next_pos
            .map(|p| i32::from(p) - 1)
            .unwrap_or(UNMAPPED_POSITION);

        let read_name = if self.read_name.is_empty() {
            b"*\x00".to_vec()
        } else {
            self.read_name
        };

        Record {
            ref_id,
            pos,
            mapq: self.mapq,
            bin: 0, // FIXME
            flag: self.flag,
            next_ref_id,
            next_pos,
            tlen: self.tlen,
            read_name,
            cigar: self.cigar,
            seq: self.seq,
            qual: self.qual,
            data: self.data,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        use sam::record::{Flags, MappingQuality};

        Self {
            ref_id: None,
            pos: None,
            mapq: MappingQuality::default(),
            flag: Flags::UNMAPPED,
            next_ref_id: None,
            next_pos: None,
            tlen: 0,
            read_name: Vec::new(),
            cigar: Cigar::default(),
            seq: Sequence::default(),
            qual: QualityScores::default(),
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

        assert!(builder.ref_id.is_none());
        assert!(builder.pos.is_none());
        assert!(builder.mapq.is_none());
        assert_eq!(builder.flag, sam::record::Flags::UNMAPPED);
        assert!(builder.next_ref_id.is_none());
        assert!(builder.next_pos.is_none());
        assert_eq!(builder.tlen, 0);
        assert!(builder.read_name.is_empty());
        assert!(builder.cigar.is_empty());
        assert!(builder.seq.is_empty());
        assert!(builder.qual.is_empty());
        assert!(builder.data.is_empty());
    }
}
