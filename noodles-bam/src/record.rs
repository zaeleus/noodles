//! BAM record.

mod bounds;
mod cigar;
pub mod codec;
pub mod data;
mod flags;
mod mapping_quality;
mod name;
mod position;
mod quality_scores;
mod reference_sequence_id;
mod sequence;
mod template_length;

use std::{fmt, io, mem};

use noodles_core as core;
use noodles_sam as sam;

use self::bounds::Bounds;
pub use self::{
    cigar::Cigar, data::Data, flags::Flags, mapping_quality::MappingQuality, name::Name,
    position::Position, quality_scores::QualityScores, reference_sequence_id::ReferenceSequenceId,
    sequence::Sequence, template_length::TemplateLength,
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
    pub fn reference_sequence_id(&self) -> Option<ReferenceSequenceId> {
        let src = &self.buf[bounds::REFERENCE_SEQUENCE_ID_RANGE];
        // SAFETY: `src` is 4 bytes.
        get_reference_sequence_id(src.try_into().unwrap())
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
    pub fn alignment_start(&self) -> Option<Position> {
        let src = &self.buf[bounds::ALIGNMENT_START_RANGE];
        // SAFETY: `src` is 4 bytes.
        get_position(src.try_into().unwrap())
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
        const MISSING: u8 = 255;

        match self.buf[bounds::MAPPING_QUALITY_INDEX] {
            MISSING => None,
            n => Some(MappingQuality::new(n)),
        }
    }

    /// Returns the flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::record::Flags;
    /// let record = bam::Record::default();
    /// assert_eq!(Flags::from(record.flags()), Flags::UNMAPPED);
    /// ```
    pub fn flags(&self) -> Flags {
        let src = &self.buf[bounds::FLAGS_RANGE];
        // SAFETY: `src` is 2 bytes.
        let n = u16::from_le_bytes(src.try_into().unwrap());
        Flags::new(n)
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
    pub fn mate_reference_sequence_id(&self) -> Option<ReferenceSequenceId> {
        let src = &self.buf[bounds::MATE_REFERENCE_SEQUENCE_ID_RANGE];
        // SAFETY: `src` is 4 bytes.
        get_reference_sequence_id(src.try_into().unwrap())
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
    pub fn mate_alignment_start(&self) -> Option<Position> {
        let src = &self.buf[bounds::MATE_ALIGNMENT_START_RANGE];
        get_position(src.try_into().unwrap())
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
    pub fn template_length(&self) -> TemplateLength {
        let src = &self.buf[bounds::TEMPLATE_LENGTH_RANGE];
        // SAFETY: `src` is 4 bytes.
        let n = i32::from_le_bytes(src.try_into().unwrap());
        TemplateLength::new(n)
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
        const MISSING: &[u8] = &[b'*', 0x00];

        match &self.buf[self.bounds.name_range()] {
            MISSING => None,
            buf => Some(Name::new(buf)),
        }
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
        let src = &self.buf[self.bounds.cigar_range()];
        Cigar::new(src)
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
        let src = &self.buf[self.bounds.sequence_range()];
        let quality_scores_range = self.bounds.quality_scores_range();
        let base_count = quality_scores_range.end - quality_scores_range.start;
        Sequence::new(src, base_count)
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
        let src = &self.buf[self.bounds.quality_scores_range()];
        QualityScores::new(src)
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
        let src = &self.buf[self.bounds.data_range()];
        Data::new(src)
    }

    pub(crate) fn index(&mut self) -> io::Result<()> {
        index(&self.buf[..], &mut self.bounds)
    }
}

fn get_reference_sequence_id(src: [u8; 4]) -> Option<ReferenceSequenceId> {
    const UNMAPPED: i32 = -1;

    match i32::from_le_bytes(src) {
        UNMAPPED => None,
        n => Some(ReferenceSequenceId::new(n)),
    }
}

fn get_position(src: [u8; 4]) -> Option<Position> {
    const MISSING: i32 = -1;

    match i32::from_le_bytes(src) {
        MISSING => None,
        n => Some(Position::new(n)),
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
            // .field("data", &self.data()) // TODO
            .finish()
    }
}

impl sam::alignment::Record for Record {
    fn name(&self) -> Option<Box<dyn sam::alignment::record::Name + '_>> {
        let name = self.name()?;
        Some(Box::new(name))
    }

    fn flags(&self) -> Box<dyn sam::alignment::record::Flags + '_> {
        Box::new(self.flags())
    }

    fn reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        _: &'h sam::Header,
    ) -> Option<Box<dyn sam::alignment::record::ReferenceSequenceId + 'r>> {
        let reference_sequence_id = self.reference_sequence_id()?;
        Some(Box::new(reference_sequence_id))
    }

    fn alignment_start(&self) -> Option<Box<dyn sam::alignment::record::Position + '_>> {
        let alignment_start = self.alignment_start()?;
        Some(Box::new(alignment_start))
    }

    fn mapping_quality(&self) -> Option<Box<dyn sam::alignment::record::MappingQuality + '_>> {
        let mapping_quality = self.mapping_quality()?;
        Some(Box::new(mapping_quality))
    }

    fn cigar(&self, header: &sam::Header) -> Box<dyn sam::alignment::record::Cigar + '_> {
        use bytes::Buf;

        use self::data::get_raw_cigar;

        const SKIP: u8 = 3;
        const SOFT_CLIP: u8 = 4;

        fn decode_op(n: u32) -> (u8, usize) {
            ((n & 0x0f) as u8, usize::try_from(n >> 4).unwrap())
        }

        let mut src = &self.buf[self.bounds.cigar_range()];

        if src.len() == 2 * mem::size_of::<u32>() {
            let k = self.sequence().len();

            let Some(Ok(m)) = self
                .reference_sequence(header)
                .map(|result| result.map(|(_, rs)| rs.length().get()))
            else {
                return Box::new(self.cigar());
            };

            // SAFETY: `src` is 8 bytes.
            let op_1 = decode_op(src.get_u32_le());
            let op_2 = decode_op(src.get_u32_le());

            if op_1 == (SOFT_CLIP, k) && op_2 == (SKIP, m) {
                let mut data_src = &self.buf[self.bounds.data_range()];

                if let Ok(Some(buf)) = get_raw_cigar(&mut data_src) {
                    return Box::new(Cigar::new(buf));
                }
            }
        }

        Box::new(self.cigar())
    }

    fn mate_reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        _: &'h sam::Header,
    ) -> Option<Box<dyn sam::alignment::record::ReferenceSequenceId + 'r>> {
        let mate_reference_sequence_id = self.mate_reference_sequence_id()?;
        Some(Box::new(mate_reference_sequence_id))
    }

    fn mate_alignment_start(&self) -> Option<Box<dyn sam::alignment::record::Position + '_>> {
        let mate_alignment_start = self.mate_alignment_start()?;
        Some(Box::new(mate_alignment_start))
    }

    fn template_length(&self) -> Box<dyn sam::alignment::record::TemplateLength + '_> {
        Box::new(self.template_length())
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

        let flags = sam::record::Flags::from(lazy_record.flags());
        builder = builder.set_flags(flags);

        if let Some(reference_sequence_id) = lazy_record.reference_sequence_id() {
            let id = usize::try_from(reference_sequence_id)?;
            builder = builder.set_reference_sequence_id(id);
        }

        if let Some(alignment_start) = lazy_record.alignment_start() {
            let start = core::Position::try_from(alignment_start)?;
            builder = builder.set_alignment_start(start);
        }

        if let Some(mapping_quality) = lazy_record.mapping_quality() {
            // SAFETY: `mapping_quality` cannot be missing.
            let mapping_quality =
                sam::record::MappingQuality::new(u8::from(mapping_quality)).unwrap();
            builder = builder.set_mapping_quality(mapping_quality);
        }

        builder = builder.set_cigar(lazy_record.cigar().try_into()?);

        if let Some(mate_reference_sequence_id) = lazy_record.mate_reference_sequence_id() {
            let id = usize::try_from(mate_reference_sequence_id)?;
            builder = builder.set_mate_reference_sequence_id(id);
        }

        if let Some(mate_alignment_start) = lazy_record.mate_alignment_start() {
            let start = core::Position::try_from(mate_alignment_start)?;
            builder = builder.set_mate_alignment_start(start);
        }

        let template_length = i32::from(lazy_record.template_length());
        builder = builder.set_template_length(template_length);

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

        use noodles_core::Position;
        use noodles_sam::{
            alignment::{
                record::{
                    cigar::{op::Kind, Op},
                    data::field::tag,
                },
                record_buf::{data::field::Value, Cigar, Sequence},
                RecordBuf,
            },
            header::record::value::{map::ReferenceSequence, Map},
            record::Flags,
        };

        use crate::record::codec::encode;

        const BASE_COUNT: usize = 65536;

        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(131072) {
            Some(n) => n,
            None => unreachable!(),
        };

        let mut buf = Vec::new();

        let header = sam::Header::builder()
            .add_reference_sequence("sq0".parse()?, Map::<ReferenceSequence>::new(SQ0_LN))
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
                [(tag::ALIGNMENT_HIT_COUNT, Value::from(1))]
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

        let cigar = sam::alignment::Record::cigar(&record, &header);
        assert_eq!(cigar.len(), BASE_COUNT);

        Ok(())
    }
}
