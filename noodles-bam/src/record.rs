//! BAM record.

mod cigar;
pub mod codec;
pub mod data;
mod quality_scores;
mod sequence;

use std::{fmt, io};

use bstr::BStr;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::record::{Flags, MappingQuality},
};

pub use self::{cigar::Cigar, data::Data, quality_scores::QualityScores, sequence::Sequence};
use super::RecordRef;

/// A BAM record.
#[derive(Clone, Eq, PartialEq)]
pub struct Record(pub(crate) Vec<u8>);

impl Record {
    /// Returns the reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.reference_sequence_id().is_none());
    /// ```
    pub fn reference_sequence_id(&self) -> Option<io::Result<usize>> {
        self.as_record_ref().reference_sequence_id()
    }

    /// Returns the alignment start.
    ///
    /// This position is 1-based, inclusive.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.alignment_start().is_none());
    /// ```
    pub fn alignment_start(&self) -> Option<io::Result<Position>> {
        self.as_record_ref().alignment_start()
    }

    /// Returns the mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.mapping_quality().is_none());
    /// ```
    pub fn mapping_quality(&self) -> Option<MappingQuality> {
        self.as_record_ref().mapping_quality()
    }

    /// Returns the flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::alignment::record::Flags;
    /// let record = bam::Record::default();
    /// assert_eq!(Flags::from(record.flags()), Flags::UNMAPPED);
    /// ```
    pub fn flags(&self) -> Flags {
        self.as_record_ref().flags()
    }

    /// Returns the mate reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.mate_reference_sequence_id().is_none());
    /// ```
    pub fn mate_reference_sequence_id(&self) -> Option<io::Result<usize>> {
        self.as_record_ref().mate_reference_sequence_id()
    }

    /// Returns the mate alignment start.
    ///
    /// This position is 1-based, inclusive.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.mate_alignment_start().is_none());
    /// ```
    pub fn mate_alignment_start(&self) -> Option<io::Result<Position>> {
        self.as_record_ref().mate_alignment_start()
    }

    /// Returns the template length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert_eq!(i32::from(record.template_length()), 0);
    /// ```
    pub fn template_length(&self) -> i32 {
        self.as_record_ref().template_length()
    }

    /// Returns the read name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.name().is_none());
    /// ```
    pub fn name(&self) -> Option<&BStr> {
        self.as_record_ref().name()
    }

    /// Returns the CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.cigar().is_empty());
    /// ```
    pub fn cigar(&self) -> Cigar<'_> {
        self.as_record_ref().cigar()
    }

    /// Returns the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.sequence().is_empty());
    /// ```
    pub fn sequence(&self) -> Sequence<'_> {
        self.as_record_ref().sequence()
    }

    /// Returns the quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.quality_scores().is_empty());
    /// ```
    pub fn quality_scores(&self) -> QualityScores<'_> {
        self.as_record_ref().quality_scores()
    }

    /// Returns the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.data().is_empty());
    /// ```
    pub fn data(&self) -> Data<'_> {
        self.as_record_ref().data()
    }

    fn as_record_ref(&self) -> RecordRef<'_> {
        // SAFETY: The inner buffer was validated to look record-like on read.
        RecordRef::new_unchecked(&self.0)
    }
}

impl Default for Record {
    fn default() -> Self {
        Self(vec![
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x02, // l_read_name = 2
            0xff, // mapq = 255
            0x48, 0x12, // bin = 4680
            0x00, 0x00, // n_cigar_op = 0
            0x04, 0x00, // flag = 4
            0x00, 0x00, 0x00, 0x00, // l_seq = 0
            0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
            0xff, 0xff, 0xff, 0xff, // next_pos = -1
            0x00, 0x00, 0x00, 0x00, // tlen = 0
            b'*', 0x00, // read_name = "*\x00"
        ])
    }
}

impl fmt::Debug for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Record")
            .field("reference_sequence_id", &self.reference_sequence_id())
            .field("alignment_start", &self.alignment_start())
            .field("mapping_quality", &self.mapping_quality())
            .field("flags", &self.flags())
            .field(
                "mate_reference_sequence_id",
                &self.mate_reference_sequence_id(),
            )
            .field("mate_alignment_start", &self.mate_alignment_start())
            .field("template_length", &self.template_length())
            .field("name", &self.name())
            .field("cigar", &self.cigar())
            .field("sequence", &self.sequence())
            .field("quality_scores", &self.quality_scores())
            .field("data", &self.data())
            .finish()
    }
}

impl sam::alignment::Record for Record {
    fn name(&self) -> Option<&BStr> {
        self.name()
    }

    fn flags(&self) -> io::Result<sam::alignment::record::Flags> {
        Ok(self.flags())
    }

    fn reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        _: &'h sam::Header,
    ) -> Option<io::Result<usize>> {
        self.reference_sequence_id()
    }

    fn alignment_start(&self) -> Option<io::Result<Position>> {
        self.alignment_start()
    }

    fn mapping_quality(&self) -> Option<io::Result<sam::alignment::record::MappingQuality>> {
        self.mapping_quality().map(Ok)
    }

    fn cigar(&self) -> Box<dyn sam::alignment::record::Cigar + '_> {
        Box::new(self.cigar())
    }

    fn cigar_ref(&self) -> sam::alignment::record::CigarRef<'_> {
        let src = self.cigar().as_bytes();
        sam::alignment::record::CigarRef::FourBytePacked(src)
    }

    fn mate_reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        _: &'h sam::Header,
    ) -> Option<io::Result<usize>> {
        self.mate_reference_sequence_id()
    }

    fn mate_alignment_start(&self) -> Option<io::Result<Position>> {
        self.mate_alignment_start()
    }

    fn template_length(&self) -> io::Result<i32> {
        Ok(self.template_length())
    }

    fn sequence(&self) -> Box<dyn sam::alignment::record::Sequence + '_> {
        Box::new(self.sequence())
    }

    fn sequence_ref(&self) -> sam::alignment::record::SequenceRef<'_> {
        use sam::alignment::record::sequence_ref::FourBitPacked;

        let sequence = self.sequence();
        let src = sequence.as_bytes();
        let base_count = sequence.len();
        sam::alignment::record::SequenceRef::FourBitPacked(FourBitPacked::new(src, base_count))
    }

    fn quality_scores(&self) -> Box<dyn sam::alignment::record::QualityScores + '_> {
        Box::new(self.quality_scores())
    }

    fn quality_scores_ref(&self) -> sam::alignment::record::QualityScoresRef<'_> {
        let src = self.quality_scores().as_bytes();
        sam::alignment::record::QualityScoresRef::Raw(src)
    }

    fn data(&self) -> Box<dyn sam::alignment::record::Data<'_> + '_> {
        Box::new(self.data())
    }

    fn data_ref(&self) -> sam::alignment::record::DataRef<'_> {
        let src = self.data().as_bytes();
        sam::alignment::record::DataRef::FieldEncoded(src)
    }
}

pub(super) fn try_to_reference_sequence_id(n: i32) -> io::Result<usize> {
    usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

pub(super) fn try_to_position(n: i32) -> io::Result<Position> {
    usize::try_from(n)
        .map(|m| m + 1)
        .and_then(Position::try_from)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_with_oversized_cigar() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZero;

        use noodles_sam::{
            alignment::{
                RecordBuf,
                record::{
                    Flags,
                    cigar::{Op, op::Kind},
                    data::field::Tag,
                },
                record_buf::{Cigar as CigarBuf, Sequence, data::field::Value},
            },
            header::record::value::{Map, map::ReferenceSequence},
        };

        use crate::record::codec::encode;

        const BASE_COUNT: usize = 65536;

        let mut buf = Vec::new();

        let header = sam::Header::builder()
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(const { NonZero::new(131072).unwrap() }),
            )
            .build();

        let cigar = CigarBuf::from(vec![Op::new(Kind::Match, 1); BASE_COUNT]);
        let sequence = Sequence::from(vec![b'A'; BASE_COUNT]);

        let record = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::MIN)
            .set_cigar(cigar)
            .set_sequence(sequence)
            .set_data(
                [(Tag::ALIGNMENT_HIT_COUNT, Value::from(1))]
                    .into_iter()
                    .collect(),
            )
            .build();

        encode(&mut buf, &header, &record)?;

        let record = Record(buf);
        assert_eq!(record.cigar().len(), BASE_COUNT);

        Ok(())
    }
}
