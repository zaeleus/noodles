use std::{
    io, mem,
    num::NonZeroUsize,
    ops::{Range, RangeFrom},
};

use byteorder::{ByteOrder, LittleEndian};
use bytes::Buf;
use noodles_core::Position;
use noodles_sam as sam;

const REFERENCE_SEQUENCE_ID_RANGE: Range<usize> = 0..4;
const ALIGNMENT_START_RANGE: Range<usize> = 4..8;
const MAPPING_QUALITY_RANGE: Range<usize> = 9..10;
const FLAGS_RANGE: Range<usize> = 14..16;
const MATE_REFERENCE_SEQUENCE_ID_RANGE: Range<usize> = 20..24;
const MATE_ALIGNMENT_START_RANGE: Range<usize> = 24..28;
const TEMPLATE_LENGTH_RANGE: Range<usize> = 28..32;

#[derive(Clone, Debug, Eq, PartialEq)]
struct Bounds {
    read_name_end: usize,
    cigar_end: usize,
    sequence_end: usize,
    quality_scores_end: usize,
}

impl Bounds {
    fn read_name_range(&self) -> Range<usize> {
        TEMPLATE_LENGTH_RANGE.end..self.read_name_end
    }

    fn cigar_range(&self) -> Range<usize> {
        self.read_name_end..self.cigar_end
    }

    fn sequence_range(&self) -> Range<usize> {
        self.cigar_end..self.sequence_end
    }

    fn quality_scores_range(&self) -> Range<usize> {
        self.sequence_end..self.quality_scores_end
    }

    fn data_range(&self) -> RangeFrom<usize> {
        self.quality_scores_end..
    }
}

/// An immutable, lazily-evalulated BAM record.
///
/// The fields are _not_ memoized.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LazyRecord {
    pub(crate) buf: Vec<u8>,
    bounds: Bounds,
}

impl LazyRecord {
    /// Returns the reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::LazyRecord::default();
    /// assert!(record.reference_sequence_id()?.is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn reference_sequence_id(&self) -> io::Result<Option<usize>> {
        use crate::reader::record::get_reference_sequence_id;
        let mut src = &self.buf[REFERENCE_SEQUENCE_ID_RANGE];
        get_reference_sequence_id(&mut src)
    }

    /// Returns the alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::LazyRecord::default();
    /// assert!(record.alignment_start()?.is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn alignment_start(&self) -> io::Result<Option<Position>> {
        use crate::reader::record::get_position;
        let mut src = &self.buf[ALIGNMENT_START_RANGE];
        get_position(&mut src)
    }

    /// Returns the mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::LazyRecord::default();
    /// assert!(record.mapping_quality()?.is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn mapping_quality(&self) -> io::Result<Option<sam::record::MappingQuality>> {
        use crate::reader::record::get_mapping_quality;
        let mut src = &self.buf[MAPPING_QUALITY_RANGE];
        get_mapping_quality(&mut src)
    }

    /// Returns the flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::record::Flags;
    ///
    /// let record = bam::LazyRecord::default();
    /// assert_eq!(record.flags()?, Flags::UNMAPPED);
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn flags(&self) -> io::Result<sam::record::Flags> {
        use crate::reader::record::get_flags;
        let mut src = &self.buf[FLAGS_RANGE];
        get_flags(&mut src)
    }

    /// Returns the mate reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::LazyRecord::default();
    /// assert!(record.mate_reference_sequence_id()?.is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn mate_reference_sequence_id(&self) -> io::Result<Option<usize>> {
        use crate::reader::record::get_reference_sequence_id;
        let mut src = &self.buf[MATE_REFERENCE_SEQUENCE_ID_RANGE];
        get_reference_sequence_id(&mut src)
    }

    /// Returns the mate alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::LazyRecord::default();
    /// assert!(record.mate_alignment_start()?.is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn mate_alignment_start(&self) -> io::Result<Option<Position>> {
        use crate::reader::record::get_position;
        let mut src = &self.buf[MATE_ALIGNMENT_START_RANGE];
        get_position(&mut src)
    }

    /// Returns the template length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::LazyRecord::default();
    /// assert_eq!(record.template_length(), 0);
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn template_length(&self) -> i32 {
        let src = &self.buf[TEMPLATE_LENGTH_RANGE];
        LittleEndian::read_i32(src)
    }

    /// Returns the read name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::LazyRecord::default();
    /// assert!(record.read_name()?.is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_name(&self) -> io::Result<Option<sam::record::ReadName>> {
        use crate::reader::record::get_read_name;

        let mut src = &self.buf[self.bounds.read_name_range()];
        let mut read_name = None;
        let len = NonZeroUsize::try_from(src.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        get_read_name(&mut src, &mut read_name, len)?;

        Ok(read_name)
    }

    /// Returns the CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::LazyRecord::default();
    /// assert!(record.cigar()?.is_empty());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn cigar(&self) -> io::Result<sam::record::Cigar> {
        use crate::reader::record::get_cigar;

        let mut src = &self.buf[self.bounds.cigar_range()];
        let mut cigar = sam::record::Cigar::default();
        let op_count = src.len() / 4;
        get_cigar(&mut src, &mut cigar, op_count)?;

        Ok(cigar)
    }

    /// Returns the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::LazyRecord::default();
    /// assert!(record.sequence()?.is_empty());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn sequence(&self) -> io::Result<sam::record::Sequence> {
        use crate::reader::record::get_sequence;

        let mut src = &self.buf[self.bounds.sequence_range()];
        let mut sequence = sam::record::Sequence::default();

        let quality_scores_range = self.bounds.quality_scores_range();
        let base_count = quality_scores_range.end - quality_scores_range.start;

        get_sequence(&mut src, &mut sequence, base_count)?;

        Ok(sequence)
    }

    /// Returns the quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::LazyRecord::default();
    /// assert!(record.quality_scores()?.is_empty());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn quality_scores(&self) -> io::Result<sam::record::QualityScores> {
        use crate::reader::record::get_quality_scores;

        let mut src = &self.buf[self.bounds.quality_scores_range()];
        let mut quality_scores = sam::record::QualityScores::default();
        let base_count = src.len();
        get_quality_scores(&mut src, &mut quality_scores, base_count)?;

        Ok(quality_scores)
    }

    /// Returns the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::LazyRecord::default();
    /// assert!(record.data()?.is_empty());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn data(&self) -> io::Result<sam::record::Data> {
        use crate::reader::record::get_data;

        let mut src = &self.buf[self.bounds.data_range()];
        let mut data = sam::record::Data::default();
        get_data(&mut src, &mut data)?;

        Ok(data)
    }

    pub(crate) fn index(&mut self) -> io::Result<()> {
        index(&self.buf[..], &mut self.bounds)
    }
}

impl Default for LazyRecord {
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
            read_name_end: buf.len(),
            cigar_end: buf.len(),
            sequence_end: buf.len(),
            quality_scores_end: buf.len(),
        };

        Self { buf, bounds }
    }
}

fn index(buf: &[u8], bounds: &mut Bounds) -> io::Result<()> {
    const MIN_BUF_LENGTH: usize = TEMPLATE_LENGTH_RANGE.end;
    const READ_NAME_LENGTH_RANGE: Range<usize> = 8..9;
    const CIGAR_OP_COUNT_RANGE: Range<usize> = 12..14;
    const READ_LENGTH_RANGE: Range<usize> = 16..20;

    if buf.len() < MIN_BUF_LENGTH {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let mut src = &buf[READ_NAME_LENGTH_RANGE];
    let l_read_name = NonZeroUsize::new(usize::from(src.get_u8()))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid l_read_name"))?;

    let mut src = &buf[CIGAR_OP_COUNT_RANGE];
    let n_cigar_op = usize::from(src.get_u16_le());

    let mut src = &buf[READ_LENGTH_RANGE];
    let l_seq = usize::try_from(src.get_u32_le())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let mut i = TEMPLATE_LENGTH_RANGE.end;
    i += usize::from(l_read_name);
    bounds.read_name_end = i;

    i += mem::size_of::<u32>() * n_cigar_op;
    bounds.cigar_end = i;

    i += (l_seq + 1) / 2;
    bounds.sequence_end = i;

    i += l_seq;
    bounds.quality_scores_end = i;

    if buf.len() < i {
        Err(io::Error::from(io::ErrorKind::UnexpectedEof))
    } else {
        Ok(())
    }
}
