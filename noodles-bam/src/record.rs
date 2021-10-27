//! BAM record and fields.

pub mod cigar;
mod convert;
pub mod data;
pub mod quality_scores;
pub mod reference_sequence_id;
pub mod sequence;

pub use self::{
    cigar::Cigar, data::Data, quality_scores::QualityScores,
    reference_sequence_id::ReferenceSequenceId, sequence::Sequence,
};

use std::{
    ffi::{self, CStr},
    fmt, mem,
};

use noodles_sam as sam;

pub(crate) const UNMAPPED_POSITION: i32 = -1;

/// A BAM record.
///
/// A BAM record encodes the same fields as a SAM record:
///
///  * reference sequence ID (`RNAME` equiv.),
///  * position (`POS`),
///  * mapping quality (`MAPQ`),
///  * flags (`FLAG`),
///  * mate reference sequence ID (`RNEXT` equiv.),
///  * mate position (`PNEXT`),
///  * template length (`TLEN`),
///  * read name (`QNAME`),
///  * CIGAR operations (`CIGAR`),
///  * sequence (`SEQ`),
///  * quality scores (`QUAL`), and
///  * optional data fields.
///
/// Additionally, it encodes the BAM index bin (`bin`).
///
/// A `bam::Record` stores raw values, and the fields should be considered immutable.
#[derive(Clone, Eq, PartialEq)]
pub struct Record {
    pub(crate) ref_id: i32,
    pub(crate) pos: i32,
    pub(crate) mapq: u8,
    pub(crate) bin: u16,
    pub(crate) flag: u16,
    pub(crate) l_seq: usize,
    pub(crate) next_ref_id: i32,
    pub(crate) next_pos: i32,
    pub(crate) tlen: i32,
    pub(crate) read_name: Vec<u8>,
    cigar: Cigar,
    pub(crate) seq: Vec<u8>,
    pub(crate) qual: Vec<u8>,
    pub(crate) data: Vec<u8>,
}

impl Record {
    pub(crate) fn block_size(&self) -> usize {
        mem::size_of::<i32>() // ref_id
            + mem::size_of::<i32>() // pos
            + mem::size_of::<u8>() // l_read_name
            + mem::size_of::<u8>() // mapq
            + mem::size_of::<u16>() // bin
            + mem::size_of::<u16>() // n_cigar_op
            + mem::size_of::<u16>() // flag
            + mem::size_of::<u32>() // l_seq
            + mem::size_of::<i32>() // next_ref_id
            + mem::size_of::<i32>() // next_pos
            + mem::size_of::<i32>() // tlen
            + self.read_name.len()
            + (mem::size_of::<u32>() * self.cigar.len())
            + self.seq.len()
            + self.qual.len()
            + self.data.len()
    }

    /// Returns the reference sequence ID of this record.
    ///
    /// The reference sequence ID is the index of the associated reference sequence in the SAM
    /// header or BAM reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.reference_sequence_id().is_none());
    /// ```
    pub fn reference_sequence_id(&self) -> Option<ReferenceSequenceId> {
        let id = self.ref_id;

        if id == reference_sequence_id::UNMAPPED {
            None
        } else {
            ReferenceSequenceId::try_from(id).ok()
        }
    }

    /// Returns the start position of this record.
    ///
    /// Despite the BAM format using 0-based positions, this normalizes the value as a 1-based
    /// position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.position().is_none());
    /// ```
    pub fn position(&self) -> Option<sam::record::Position> {
        let pos = self.pos;

        if pos == UNMAPPED_POSITION {
            None
        } else {
            sam::record::Position::try_from(pos + 1).ok()
        }
    }

    /// Returns the mapping quality of this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam as sam;
    /// let record = bam::Record::default();
    /// assert!(record.mapping_quality().is_none());
    /// ```
    pub fn mapping_quality(&self) -> sam::record::MappingQuality {
        sam::record::MappingQuality::from(self.mapq)
    }

    /// Returns the index bin that includes this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert_eq!(record.bin(), 4680);
    /// ```
    pub fn bin(&self) -> u16 {
        self.bin
    }

    /// Returns the SAM flags of this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam as sam;
    /// let record = bam::Record::default();
    /// assert_eq!(record.flags(), sam::record::Flags::UNMAPPED);
    /// ```
    pub fn flags(&self) -> sam::record::Flags {
        sam::record::Flags::from(self.flag)
    }

    /// Returns the reference sequence ID of the mate of this record.
    ///
    /// The mate reference sequence ID is the index of the associated reference sequence in the SAM
    /// header or BAM reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.mate_reference_sequence_id().is_none());
    /// ```
    pub fn mate_reference_sequence_id(&self) -> Option<ReferenceSequenceId> {
        let id = self.next_ref_id;

        if id == reference_sequence_id::UNMAPPED {
            None
        } else {
            ReferenceSequenceId::try_from(id).ok()
        }
    }

    /// Returns the start position of the mate of this record.
    ///
    /// Despite the BAM format using 0-based positions, this normalizes the value as a 1-based
    /// position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.mate_position().is_none());
    /// ```
    pub fn mate_position(&self) -> Option<sam::record::Position> {
        let pos = self.next_pos;

        if pos == UNMAPPED_POSITION {
            None
        } else {
            sam::record::Position::try_from(pos + 1).ok()
        }
    }

    /// Returns the template length of this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert_eq!(record.template_length(), 0);
    /// ```
    pub fn template_length(&self) -> i32 {
        self.tlen
    }

    /// Returns the read name of this record.
    ///
    /// This is also called the query name.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::ffi;
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert_eq!(record.read_name()?.to_bytes(), b"*");
    /// # Ok::<(), ffi::FromBytesWithNulError>(())
    /// ```
    pub fn read_name(&self) -> Result<&CStr, ffi::FromBytesWithNulError> {
        CStr::from_bytes_with_nul(&self.read_name)
    }

    /// Returns the CIGAR operations that describe how the read was mapped.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.cigar().is_empty());
    /// ```
    pub fn cigar(&self) -> &Cigar {
        &self.cigar
    }

    /// Returns a mutable reference to the CIGAR.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::{self as bam, record::cigar::Op};
    /// use noodles_sam::record::cigar::op::Kind;
    ///
    /// let mut record = bam::Record::default();
    ///
    /// let raw_op = Op::new(Kind::Match, 36).map(u32::from)?;
    /// record.cigar_mut().push(raw_op);
    ///
    /// assert_eq!(**record.cigar(), [raw_op]);
    /// Ok::<_, bam::record::cigar::op::LengthError>(())
    /// ```
    pub fn cigar_mut(&mut self) -> &mut Cigar {
        &mut self.cigar
    }

    /// Returns the bases in the sequence of this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.sequence().is_empty());
    /// ```
    pub fn sequence(&self) -> Sequence<'_> {
        Sequence::new(&self.seq, self.l_seq)
    }

    /// Returns the quality score for each base in the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.quality_scores().is_empty());
    /// ```
    pub fn quality_scores(&self) -> QualityScores<'_> {
        QualityScores::new(&self.qual)
    }

    /// Returns the optional data fields for this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert!(record.data().is_empty());
    /// ```
    pub fn data(&self) -> Data<'_> {
        Data::new(&self.data)
    }
}

impl Default for Record {
    fn default() -> Self {
        Self {
            ref_id: -1,
            pos: -1,
            mapq: 255,
            bin: 4680,
            flag: 0x0004,
            l_seq: 0,
            next_ref_id: -1,
            next_pos: -1,
            tlen: 0,
            read_name: b"*\x00".to_vec(),
            cigar: Cigar::default(),
            seq: Vec::new(),
            qual: Vec::new(),
            data: Vec::new(),
        }
    }
}

impl fmt::Debug for Record {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.debug_struct("Record")
            .field("block_size", &self.block_size())
            .field("ref_id", &self.reference_sequence_id())
            .field("pos", &self.position())
            .field("mapq", &self.mapping_quality())
            .field("bin", &self.bin())
            .field("flag", &self.flags())
            .field("l_seq", &self.l_seq)
            .field("next_ref_id", &self.mate_reference_sequence_id())
            .field("next_pos", &self.mate_position())
            .field("tlen", &self.template_length())
            .field("read_name", &self.read_name())
            .field("cigar", &self.cigar())
            .field("seq", &self.sequence())
            .field("data", &self.data())
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_record() -> Record {
        use sam::record::Flags;

        Record {
            ref_id: 10,
            pos: 61061,
            mapq: 12,
            bin: 4684,
            flag: u16::from(Flags::PAIRED | Flags::READ_1),
            l_seq: 4,
            next_ref_id: 10,
            next_pos: 61152,
            tlen: 166,
            read_name: b"r0\x00".to_vec(),
            cigar: Cigar::from(vec![0x00000040]), // 4M
            seq: vec![0x18, 0x42],                // ATGC
            qual: vec![0x1f, 0x1d, 0x1e, 0x20],   // @>?A
            data: vec![
                0x4e, 0x4d, 0x43, 0x00, // NM:i:0
                0x50, 0x47, 0x5a, 0x53, 0x4e, 0x41, 0x50, 0x00, // PG:Z:SNAP
            ],
        }
    }

    #[test]
    fn test_block_size() {
        let record = build_record();
        assert_eq!(record.block_size(), 57);
    }

    #[test]
    fn test_reference_sequence_id() {
        let record = build_record();
        assert_eq!(record.reference_sequence_id().map(i32::from), Some(10));
    }

    #[test]
    fn test_position() {
        let record = build_record();
        assert_eq!(record.position().map(i32::from), Some(61062));
    }

    #[test]
    fn test_mapping_quality() {
        let record = build_record();
        assert_eq!(
            record.mapping_quality(),
            sam::record::MappingQuality::from(12)
        );
    }

    #[test]
    fn test_bin() {
        let record = build_record();
        assert_eq!(record.bin(), 4684);
    }

    #[test]
    fn test_flags() {
        use sam::record::Flags;
        let record = build_record();
        assert_eq!(record.flags(), Flags::PAIRED | Flags::READ_1);
    }

    #[test]
    fn test_l_seq() {
        let record = build_record();
        assert_eq!(record.l_seq, 4);
    }

    #[test]
    fn test_mate_reference_sequence_id() {
        let record = build_record();
        assert_eq!(record.mate_reference_sequence_id().map(i32::from), Some(10));
    }

    #[test]
    fn test_mate_position() {
        let record = build_record();
        assert_eq!(record.mate_position().map(i32::from), Some(61153));
    }

    #[test]
    fn test_template_length() {
        let record = build_record();
        assert_eq!(record.template_length(), 166);
    }

    #[test]
    fn test_read_name() {
        let record = build_record();

        assert_eq!(
            record.read_name().map(|name| name.to_bytes()),
            Ok("r0".as_bytes())
        );
    }

    #[test]
    fn test_cigar() {
        let record = build_record();
        assert_eq!(**record.cigar(), [0x00000040]);
    }

    #[test]
    fn test_sequence() {
        let record = build_record();
        assert_eq!(*record.sequence(), [0x18, 0x42]);
    }

    #[test]
    fn test_quality_scores() {
        let record = build_record();
        assert_eq!(*record.quality_scores(), [0x1f, 0x1d, 0x1e, 0x20]);
    }

    #[test]
    fn test_data() {
        let record = build_record();

        assert_eq!(
            *record.data(),
            [0x4e, 0x4d, 0x43, 0x00, 0x50, 0x47, 0x5a, 0x53, 0x4e, 0x41, 0x50, 0x00,]
        );
    }
}
