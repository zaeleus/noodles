//! BAM record and fields.

pub mod builder;
pub mod cigar;
mod convert;
pub mod data;
pub mod quality_scores;
pub mod reference_sequence_id;
pub mod sequence;

pub use self::{
    builder::Builder, cigar::Cigar, data::Data, quality_scores::QualityScores,
    reference_sequence_id::ReferenceSequenceId, sequence::Sequence,
};

use std::{fmt, io, mem};

use noodles_sam::{
    self as sam,
    header::{ReferenceSequence, ReferenceSequences},
};

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
/// A `bam::Record` and its fields store raw values and care should be taken when manipulating
/// them.
#[derive(Clone, Eq, PartialEq)]
pub struct Record {
    ref_id: Option<ReferenceSequenceId>,
    pub(crate) pos: i32,
    mapq: Option<sam::record::MappingQuality>,
    bin: u16,
    flag: sam::record::Flags,
    next_ref_id: Option<ReferenceSequenceId>,
    pub(crate) next_pos: i32,
    tlen: i32,
    read_name: Vec<u8>,
    cigar: Cigar,
    seq: Sequence,
    qual: QualityScores,
    data: Data,
}

impl Record {
    /// Creates a BAM record builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let builder = bam::Record::builder();
    /// let record = builder.build()?;
    /// # Ok::<(), bam::record::builder::BuildError>(())
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

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
            + self.read_name.len() + mem::size_of::<u8>() // read_name + NUL terminator
            + (mem::size_of::<u32>() * self.cigar.len())
            + self.seq.as_ref().len()
            + self.qual.len()
            + self.data.as_ref().len()
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
        self.ref_id
    }

    /// Returns a mutable reference to the reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::{self as bam, record::ReferenceSequenceId};
    ///
    /// let mut record = bam::Record::default();
    /// *record.reference_sequence_id_mut() = Some(ReferenceSequenceId::from(1));
    ///
    /// assert_eq!(
    ///     record.reference_sequence_id(),
    ///     Some(ReferenceSequenceId::from(1))
    /// );
    /// ```
    pub fn reference_sequence_id_mut(&mut self) -> &mut Option<ReferenceSequenceId> {
        &mut self.ref_id
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
    /// let record = bam::Record::default();
    /// assert!(record.mapping_quality().is_none());
    /// ```
    pub fn mapping_quality(&self) -> Option<sam::record::MappingQuality> {
        self.mapq
    }

    /// Returns a mutable reference to the mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::record::MappingQuality;
    ///
    /// let mut record = bam::Record::default();
    /// *record.mapping_quality_mut() = MappingQuality::try_from(13).map(Some)?;
    ///
    /// assert_eq!(record.mapping_quality().map(u8::from), Some(13));
    /// # Ok::<_, noodles_sam::record::mapping_quality::ParseError>(())
    /// ```
    pub fn mapping_quality_mut(&mut self) -> &mut Option<sam::record::MappingQuality> {
        &mut self.mapq
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

    pub(crate) fn bin_mut(&mut self) -> &mut u16 {
        &mut self.bin
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
        self.flag
    }

    /// Returns a mutable reference to the flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::record::Flags;
    /// let mut record = bam::Record::default();
    /// *record.flags_mut() = Flags::PAIRED | Flags::READ_1;
    /// assert_eq!(record.flags(), Flags::PAIRED | Flags::READ_1);
    /// ```
    pub fn flags_mut(&mut self) -> &mut sam::record::Flags {
        &mut self.flag
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
        self.next_ref_id
    }

    /// Returns a mutable reference to the mate reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::{self as bam, record::ReferenceSequenceId};
    ///
    /// let mut record = bam::Record::default();
    /// *record.mate_reference_sequence_id_mut() = Some(ReferenceSequenceId::from(1));
    ///
    /// assert_eq!(
    ///     record.mate_reference_sequence_id(),
    ///     Some(ReferenceSequenceId::from(1))
    /// );
    /// ```
    pub fn mate_reference_sequence_id_mut(&mut self) -> &mut Option<ReferenceSequenceId> {
        &mut self.next_ref_id
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

    pub(crate) fn template_length_mut(&mut self) -> &mut i32 {
        &mut self.tlen
    }

    /// Returns the read name of this record.
    ///
    /// This is also called the query name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let record = bam::Record::default();
    /// assert_eq!(record.read_name(), b"*");
    /// ```
    pub fn read_name(&self) -> &[u8] {
        &self.read_name
    }

    /// Returns a mutable reference to the read name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let mut record = bam::Record::default();
    /// *record.read_name_mut() = b"r1".to_vec();
    /// assert_eq!(record.read_name(), b"r1");
    /// ```
    pub fn read_name_mut(&mut self) -> &mut Vec<u8> {
        &mut self.read_name
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
    /// let op = Op::new(Kind::Match, 36)?;
    /// record.cigar_mut().push(op);
    ///
    /// assert_eq!(record.cigar().as_ref(), [0x00000240]);
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
    pub fn sequence(&self) -> &Sequence {
        &self.seq
    }

    /// Returns a mutable reference to the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::{self as bam, record::sequence::Base};
    ///
    /// let mut record = bam::Record::default();
    ///
    /// let sequence = record.sequence_mut();
    /// sequence.push(Base::A);
    /// sequence.set_len(1);
    ///
    /// assert_eq!(record.sequence().as_ref(), [0x10]); // A
    /// assert_eq!(record.sequence().len(), 1);
    /// ```
    pub fn sequence_mut(&mut self) -> &mut Sequence {
        &mut self.seq
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
    pub fn quality_scores(&self) -> &QualityScores {
        &self.qual
    }

    /// Returns a mutable reference to the quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::record::quality_scores::Score;
    ///
    /// let mut record = bam::Record::default();
    /// record.quality_scores_mut().push(Score::try_from(8)?);
    ///
    /// assert_eq!(record.quality_scores().as_ref(), [8]);
    /// # Ok::<_, noodles_sam::record::quality_scores::score::TryFromUByteError>(())
    /// ```
    pub fn quality_scores_mut(&mut self) -> &mut QualityScores {
        &mut self.qual
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
    pub fn data(&self) -> &Data {
        &self.data
    }

    /// Returns a mutable reference to the data.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::{self as bam, record::data::{field::Value, Field}};
    /// use noodles_sam::record::data::field::Tag;
    ///
    /// let mut record = bam::Record::default();
    ///
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::UInt8(1));
    /// record.data_mut().insert(nh).transpose()?;
    ///
    /// assert_eq!(record.data().len(), 1);
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn data_mut(&mut self) -> &mut Data {
        &mut self.data
    }
}

impl sam::AlignmentRecord for Record {
    /// Returns the associated reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{header::ReferenceSequences, AlignmentRecord};
    ///
    /// let record = bam::Record::default();
    /// let reference_sequences = ReferenceSequences::default();
    ///
    /// assert!(record.reference_sequence(&reference_sequences).is_none());
    /// ```
    fn reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>> {
        get_reference_sequence(reference_sequences, self.reference_sequence_id())
    }

    /// Returns the start position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::AlignmentRecord;
    /// let record = bam::Record::default();
    /// assert!(record.alignment_start().is_none());
    /// ```
    fn alignment_start(&self) -> Option<sam::record::Position> {
        self.position()
    }

    /// Calculates the alignment span over the reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam as bam;
    /// use noodles_sam::AlignmentRecord;
    /// let record = bam::Record::default();
    /// assert_eq!(record.alignment_span()?, 0);
    /// # Ok::<_, io::Error>(())
    /// ```
    fn alignment_span(&self) -> io::Result<u32> {
        self.cigar().reference_len()
    }

    /// Returns the associated reference sequence of the mate.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{header::ReferenceSequences, AlignmentRecord};
    ///
    /// let record = bam::Record::default();
    /// let reference_sequences = ReferenceSequences::default();
    ///
    /// assert!(record.mate_reference_sequence(&reference_sequences).is_none());
    /// ```
    fn mate_reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>> {
        get_reference_sequence(reference_sequences, self.mate_reference_sequence_id())
    }

    fn mate_alignment_start(&self) -> Option<sam::record::Position> {
        self.mate_position()
    }
}

fn get_reference_sequence(
    reference_sequences: &ReferenceSequences,
    reference_sequence_id: Option<ReferenceSequenceId>,
) -> Option<io::Result<&ReferenceSequence>> {
    reference_sequence_id.map(|reference_sequence_id| {
        let id = usize::from(reference_sequence_id);

        reference_sequences
            .get_index(id)
            .map(|(_, rs)| rs)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "invalid reference sequence ID")
            })
    })
}

impl Default for Record {
    fn default() -> Self {
        use sam::record::Flags;

        Self {
            ref_id: None,
            pos: UNMAPPED_POSITION,
            mapq: None,
            bin: 4680,
            flag: Flags::UNMAPPED,
            next_ref_id: None,
            next_pos: UNMAPPED_POSITION,
            tlen: 0,
            read_name: b"*".to_vec(),
            cigar: Cigar::default(),
            seq: Sequence::default(),
            qual: QualityScores::default(),
            data: Data::default(),
        }
    }
}

impl fmt::Debug for Record {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        use std::str;

        fmt.debug_struct("Record")
            .field("block_size", &self.block_size())
            .field("ref_id", &self.reference_sequence_id())
            .field("pos", &self.position())
            .field("mapq", &self.mapping_quality())
            .field("bin", &self.bin())
            .field("flag", &self.flags())
            .field("next_ref_id", &self.mate_reference_sequence_id())
            .field("next_pos", &self.mate_position())
            .field("tlen", &self.template_length())
            .field("read_name", &str::from_utf8(self.read_name()))
            .field("cigar", &self.cigar())
            .field("seq", &self.sequence())
            .field("data", &self.data())
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use std::io;

    use super::*;

    fn build_record() -> io::Result<Record> {
        use sam::record::{Flags, MappingQuality};

        let ref_id = Some(ReferenceSequenceId::from(10));

        let mapq = MappingQuality::try_from(12)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        Ok(Record {
            ref_id,
            pos: 61061,
            mapq,
            bin: 4684,
            flag: Flags::SEGMENTED | Flags::FIRST_SEGMENT,
            next_ref_id: ref_id,
            next_pos: 61152,
            tlen: 166,
            read_name: b"r0".to_vec(),
            cigar: Cigar::from(vec![0x00000040]),    // 4M
            seq: Sequence::new(vec![0x18, 0x42], 4), // ATGC
            qual: QualityScores::from(vec![0x1f, 0x1d, 0x1e, 0x20]), // @>?A
            data: Data::try_from(vec![
                0x4e, 0x4d, 0x43, 0x00, // NM:i:0
                0x50, 0x47, 0x5a, 0x53, 0x4e, 0x41, 0x50, 0x00, // PG:Z:SNAP
            ])?,
        })
    }

    #[test]
    fn test_block_size() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.block_size(), 57);
        Ok(())
    }

    #[test]
    fn test_reference_sequence_id() -> io::Result<()> {
        let record = build_record()?;

        assert_eq!(
            record.reference_sequence_id(),
            Some(ReferenceSequenceId::from(10))
        );

        Ok(())
    }

    #[test]
    fn test_position() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.position().map(i32::from), Some(61062));
        Ok(())
    }

    #[test]
    fn test_mapping_quality() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.mapping_quality().map(u8::from), Some(12));
        Ok(())
    }

    #[test]
    fn test_bin() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.bin(), 4684);
        Ok(())
    }

    #[test]
    fn test_flags() -> io::Result<()> {
        use sam::record::Flags;
        let record = build_record()?;
        assert_eq!(record.flags(), Flags::SEGMENTED | Flags::FIRST_SEGMENT);
        Ok(())
    }

    #[test]
    fn test_mate_reference_sequence_id() -> io::Result<()> {
        let record = build_record()?;

        assert_eq!(
            record.mate_reference_sequence_id(),
            Some(ReferenceSequenceId::from(10))
        );

        Ok(())
    }

    #[test]
    fn test_mate_position() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.mate_position().map(i32::from), Some(61153));
        Ok(())
    }

    #[test]
    fn test_template_length() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.template_length(), 166);
        Ok(())
    }

    #[test]
    fn test_read_name() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.read_name(), b"r0");
        Ok(())
    }

    #[test]
    fn test_cigar() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.cigar().as_ref(), [0x00000040]);
        Ok(())
    }

    #[test]
    fn test_sequence() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.sequence().as_ref(), [0x18, 0x42]);
        Ok(())
    }

    #[test]
    fn test_quality_scores() -> io::Result<()> {
        let record = build_record()?;
        assert_eq!(record.quality_scores().as_ref(), [0x1f, 0x1d, 0x1e, 0x20]);
        Ok(())
    }

    #[test]
    fn test_data() -> io::Result<()> {
        let record = build_record()?;

        assert_eq!(
            record.data().as_ref(),
            [0x4e, 0x4d, 0x43, 0x00, 0x50, 0x47, 0x5a, 0x53, 0x4e, 0x41, 0x50, 0x00,]
        );

        Ok(())
    }
}
