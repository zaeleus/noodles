//! BAM record.

mod bounds;
pub mod codec;
pub mod field;

use std::{fmt, io, mem};

use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::record::{Flags, MappingQuality},
};

use self::{
    bounds::Bounds,
    field::{Cigar, Data, Fields, Name, QualityScores, Sequence},
};

/// A BAM record.
#[derive(Clone, Eq, PartialEq)]
pub struct Record {
    pub(crate) buf: Vec<u8>,
    bounds: Bounds,
}

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
        self.fields()
            .reference_sequence_id()
            .map(try_to_reference_sequence_id)
    }

    /// Returns the alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.alignment_start().is_none());
    /// ```
    pub fn alignment_start(&self) -> Option<io::Result<Position>> {
        self.fields().alignment_start().map(try_to_position)
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
        self.fields()
            .mapping_quality()
            .and_then(MappingQuality::new)
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
        Flags::from(self.fields().flags())
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
        self.fields()
            .mate_reference_sequence_id()
            .map(try_to_reference_sequence_id)
    }

    /// Returns the mate alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.mate_alignment_start().is_none());
    /// ```
    pub fn mate_alignment_start(&self) -> Option<io::Result<Position>> {
        self.fields().mate_alignment_start().map(try_to_position)
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
        self.fields().template_length()
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
    pub fn name(&self) -> Option<Name> {
        self.fields().name()
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
        self.fields().cigar()
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
        self.fields().sequence()
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
        self.fields().quality_scores()
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
        self.fields().data()
    }

    fn fields(&self) -> Fields<'_> {
        Fields::new(&self.buf, &self.bounds)
    }

    pub(crate) fn index(&mut self) -> io::Result<()> {
        index(&self.buf[..], &mut self.bounds)
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
    fn name(&self) -> Option<Box<dyn sam::alignment::record::Name + '_>> {
        let name = self.name()?;
        Some(Box::new(name))
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

    fn quality_scores(&self) -> Box<dyn sam::alignment::record::QualityScores + '_> {
        Box::new(self.quality_scores())
    }

    fn data(&self) -> Box<dyn sam::alignment::record::Data + '_> {
        Box::new(self.data())
    }
}

impl AsRef<[u8]> for Record {
    fn as_ref(&self) -> &[u8] {
        &self.buf
    }
}

impl TryFrom<Vec<u8>> for Record {
    type Error = io::Error;

    fn try_from(buf: Vec<u8>) -> Result<Self, Self::Error> {
        let mut bounds = Bounds {
            name_end: 0,
            cigar_end: 0,
            sequence_end: 0,
            quality_scores_end: 0,
        };

        index(&buf, &mut bounds)?;

        Ok(Record { buf, bounds })
    }
}

impl Default for Record {
    fn default() -> Self {
        let buf = vec![
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
        ];

        let bounds = Bounds {
            name_end: buf.len(),
            cigar_end: buf.len(),
            sequence_end: buf.len(),
            quality_scores_end: buf.len(),
        };

        Self { buf, bounds }
    }
}

impl TryFrom<Record> for sam::alignment::RecordBuf {
    type Error = io::Error;

    fn try_from(lazy_record: Record) -> Result<Self, Self::Error> {
        let mut builder = Self::builder();

        if let Some(name) = lazy_record.name() {
            builder = builder.set_name(name.into());
        }

        builder = builder.set_flags(lazy_record.flags());

        if let Some(reference_sequence_id) = lazy_record.reference_sequence_id().transpose()? {
            builder = builder.set_reference_sequence_id(reference_sequence_id);
        }

        if let Some(alignment_start) = lazy_record.alignment_start().transpose()? {
            builder = builder.set_alignment_start(alignment_start);
        }

        if let Some(mapping_quality) = lazy_record.mapping_quality() {
            builder = builder.set_mapping_quality(mapping_quality);
        }

        builder = builder.set_cigar(lazy_record.cigar().try_into()?);

        if let Some(mate_reference_sequence_id) =
            lazy_record.mate_reference_sequence_id().transpose()?
        {
            builder = builder.set_mate_reference_sequence_id(mate_reference_sequence_id);
        }

        if let Some(mate_alignment_start) = lazy_record.mate_alignment_start().transpose()? {
            builder = builder.set_mate_alignment_start(mate_alignment_start);
        }

        builder = builder.set_template_length(lazy_record.template_length());

        builder = builder
            .set_sequence(lazy_record.sequence().into())
            .set_quality_scores(lazy_record.quality_scores().into())
            .set_data(lazy_record.data().try_into()?);

        Ok(builder.build())
    }
}

fn index(buf: &[u8], bounds: &mut Bounds) -> io::Result<()> {
    const MIN_BUF_LENGTH: usize = bounds::TEMPLATE_LENGTH_RANGE.end;

    if buf.len() < MIN_BUF_LENGTH {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let read_name_len = usize::from(buf[bounds::NAME_LENGTH_INDEX]);
    bounds.name_end = bounds::TEMPLATE_LENGTH_RANGE.end + read_name_len;

    let src = &buf[bounds::CIGAR_OP_COUNT_RANGE];
    // SAFETY: `src` is 2 bytes.
    let cigar_op_count = usize::from(u16::from_le_bytes(src.try_into().unwrap()));
    let cigar_len = mem::size_of::<u32>() * cigar_op_count;
    bounds.cigar_end = bounds.name_end + cigar_len;

    let src = &buf[bounds::READ_LENGTH_RANGE];
    // SAFETY: `src` is 4 bytes.
    let base_count = usize::try_from(u32::from_le_bytes(src.try_into().unwrap()))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let sequence_len = (base_count + 1) / 2;
    bounds.sequence_end = bounds.cigar_end + sequence_len;

    bounds.quality_scores_end = bounds.sequence_end + base_count;

    if buf.len() < bounds.quality_scores_end {
        Err(io::Error::from(io::ErrorKind::UnexpectedEof))
    } else {
        Ok(())
    }
}

fn try_to_reference_sequence_id(n: i32) -> io::Result<usize> {
    usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn try_to_position(n: i32) -> io::Result<Position> {
    usize::try_from(n)
        .map(|m| m + 1)
        .and_then(Position::try_from)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    static DATA: &[u8] = &[
        0xff, 0xff, 0xff, 0xff, // ref_id = -1
        0xff, 0xff, 0xff, 0xff, // pos = -1
        0x02, // l_read_name = 2
        0xff, // mapq = 255
        0x48, 0x12, // bin = 4680
        0x01, 0x00, // n_cigar_op = 1
        0x04, 0x00, // flag = 4
        0x04, 0x00, 0x00, 0x00, // l_seq = 0
        0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
        0xff, 0xff, 0xff, 0xff, // next_pos = -1
        0x00, 0x00, 0x00, 0x00, // tlen = 0
        b'*', 0x00, // read_name = "*\x00"
        0x40, 0x00, 0x00, 0x00, // cigar = 4M
        0x12, 0x48, // sequence = ACGT
        b'N', b'D', b'L', b'S', // quality scores
    ];

    #[test]
    fn test_index() -> io::Result<()> {
        let mut record = Record::default();

        record.buf.clear();
        record.buf.extend(DATA);

        record.index()?;

        assert_eq!(record.bounds.name_range(), 32..34);
        assert_eq!(record.bounds.cigar_range(), 34..38);
        assert_eq!(record.bounds.sequence_range(), 38..40);
        assert_eq!(record.bounds.quality_scores_range(), 40..44);
        assert_eq!(record.bounds.data_range(), 44..);

        Ok(())
    }

    #[test]
    fn test_try_from_vec_u8_for_record() -> io::Result<()> {
        let record = Record::try_from(DATA.to_vec())?;

        assert_eq!(record.buf, DATA);

        assert_eq!(record.bounds.name_range(), 32..34);
        assert_eq!(record.bounds.cigar_range(), 34..38);
        assert_eq!(record.bounds.sequence_range(), 38..40);
        assert_eq!(record.bounds.quality_scores_range(), 40..44);
        assert_eq!(record.bounds.data_range(), 44..);

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_sam_alignment_record() -> io::Result<()> {
        let lazy_record = Record::default();
        let actual = sam::alignment::RecordBuf::try_from(lazy_record)?;

        let expected = sam::alignment::RecordBuf::default();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_cigar_with_oversized_cigar() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZeroUsize;

        use noodles_sam::{
            alignment::{
                record::{
                    cigar::{op::Kind, Op},
                    data::field::Tag,
                    Flags,
                },
                record_buf::{data::field::Value, Cigar, Sequence},
                RecordBuf,
            },
            header::record::value::{map::ReferenceSequence, Map},
        };

        use crate::record::codec::encode;

        const BASE_COUNT: usize = 65536;

        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(131072) {
            Some(n) => n,
            None => unreachable!(),
        };

        let mut buf = Vec::new();

        let header = sam::Header::builder()
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(SQ0_LN))
            .build();

        let cigar = Cigar::from(vec![Op::new(Kind::Match, 1); BASE_COUNT]);
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

        let mut record = Record {
            buf,
            ..Default::default()
        };

        record.index()?;

        assert_eq!(record.cigar().len(), BASE_COUNT);

        Ok(())
    }
}
