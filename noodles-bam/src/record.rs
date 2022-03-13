//! BAM record and fields.

pub mod builder;
mod convert;
pub mod reference_sequence_id;

pub use self::{builder::Builder, reference_sequence_id::ReferenceSequenceId};

use std::io;

use noodles_core::Position;
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
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    ref_id: Option<ReferenceSequenceId>,
    pos: Option<Position>,
    mapq: Option<sam::record::MappingQuality>,
    flag: sam::record::Flags,
    next_ref_id: Option<ReferenceSequenceId>,
    next_pos: Option<Position>,
    tlen: i32,
    read_name: Option<sam::record::ReadName>,
    cigar: sam::record::Cigar,
    seq: sam::record::Sequence,
    qual: sam::record::QualityScores,
    data: sam::record::Data,
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
    pub fn position(&self) -> Option<Position> {
        self.pos
    }

    /// Returns a mutable reference to the start position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_core::Position;
    ///
    /// let position = Position::new(8);
    ///
    /// let mut record = bam::Record::default();
    /// *record.position_mut() = position;
    ///
    /// assert_eq!(record.position(), position);
    /// ```
    pub fn position_mut(&mut self) -> &mut Option<Position> {
        &mut self.pos
    }

    /// Returns a mutable reference to the mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{record::MappingQuality, AlignmentRecord};
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

    /// Returns a mutable reference to the flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{record::Flags, AlignmentRecord};
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
    pub fn mate_position(&self) -> Option<Position> {
        self.next_pos
    }

    /// Returns a mutable reference to the position of the mate.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_core::Position;
    ///
    /// let position = Position::new(13);
    ///
    /// let mut record = bam::Record::default();
    /// *record.mate_position_mut() = position;
    ///
    /// assert_eq!(record.mate_position(), position);
    /// ```
    pub fn mate_position_mut(&mut self) -> &mut Option<Position> {
        &mut self.next_pos
    }

    pub(crate) fn template_length_mut(&mut self) -> &mut i32 {
        &mut self.tlen
    }

    /// Returns a mutable reference to the read name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{record::ReadName, AlignmentRecord};
    ///
    /// let read_name: ReadName = "r1".parse()?;
    ///
    /// let mut record = bam::Record::default();
    /// *record.read_name_mut() = Some(read_name.clone());
    ///
    /// assert_eq!(record.read_name(), Some(&read_name));
    /// # Ok::<_, noodles_sam::record::read_name::ParseError>(())
    /// ```
    pub fn read_name_mut(&mut self) -> &mut Option<sam::record::ReadName> {
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
    pub fn cigar(&self) -> &sam::record::Cigar {
        &self.cigar
    }

    /// Returns a mutable reference to the CIGAR.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::record::Cigar;
    ///
    /// let cigar: Cigar = "36M".parse()?;
    ///
    /// let mut record = bam::Record::default();
    /// *record.cigar_mut() = cigar.clone();
    ///
    /// assert_eq!(record.cigar(), &cigar);
    /// Ok::<_, noodles_sam::record::cigar::ParseError>(())
    /// ```
    pub fn cigar_mut(&mut self) -> &mut sam::record::Cigar {
        &mut self.cigar
    }

    /// Returns a mutable reference to the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{record::Sequence, AlignmentRecord};
    ///
    /// let sequence: Sequence = "ACGT".parse()?;
    ///
    /// let mut record = bam::Record::default();
    /// *record.sequence_mut() = sequence.clone();
    ///
    /// assert_eq!(record.sequence(), &sequence);
    /// # Ok::<_, noodles_sam::record::sequence::ParseError>(())
    /// ```
    pub fn sequence_mut(&mut self) -> &mut sam::record::Sequence {
        &mut self.seq
    }

    /// Returns a mutable reference to the quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{record::{quality_scores::Score, QualityScores}, AlignmentRecord};
    ///
    /// let quality_scores: QualityScores = "NDLS".parse()?;
    ///
    /// let mut record = bam::Record::default();
    /// *record.quality_scores_mut() = quality_scores.clone();
    ///
    /// assert_eq!(record.quality_scores(), &quality_scores);
    /// # Ok::<_, noodles_sam::record::quality_scores::ParseError>(())
    /// ```
    pub fn quality_scores_mut(&mut self) -> &mut sam::record::QualityScores {
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
    pub fn data(&self) -> &sam::record::Data {
        &self.data
    }

    /// Returns a mutable reference to the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::record::data::{field::{Tag, Value}, Field};
    ///
    /// let mut record = bam::Record::default();
    ///
    /// let nh = Field::new(Tag::AlignmentHitCount, Value::UInt8(1));
    /// record.data_mut().insert(nh);
    ///
    /// assert_eq!(record.data().len(), 1);
    /// ```
    pub fn data_mut(&mut self) -> &mut sam::record::Data {
        &mut self.data
    }
}

impl sam::AlignmentRecord for Record {
    fn read_name(&self) -> Option<&sam::record::ReadName> {
        self.read_name.as_ref()
    }

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

    fn flags(&self) -> sam::record::Flags {
        self.flag
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
    fn alignment_start(&self) -> Option<Position> {
        self.position()
    }

    /// Calculates the alignment span over the reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::AlignmentRecord;
    /// let record = bam::Record::default();
    /// assert_eq!(record.alignment_span(), 0);
    /// ```
    fn alignment_span(&self) -> usize {
        self.cigar().reference_len()
    }

    fn mapping_quality(&self) -> Option<sam::record::MappingQuality> {
        self.mapq
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

    fn mate_alignment_start(&self) -> Option<Position> {
        self.mate_position()
    }

    fn template_length(&self) -> i32 {
        self.tlen
    }

    fn sequence(&self) -> &sam::record::Sequence {
        &self.seq
    }

    fn quality_scores(&self) -> &sam::record::QualityScores {
        &self.qual
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
        Self {
            ref_id: None,
            pos: None,
            mapq: None,
            flag: sam::record::Flags::UNMAPPED,
            next_ref_id: None,
            next_pos: None,
            tlen: 0,
            read_name: None,
            cigar: sam::record::Cigar::default(),
            seq: sam::record::Sequence::default(),
            qual: sam::record::QualityScores::default(),
            data: sam::record::Data::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let record = Record::default();

        assert!(record.ref_id.is_none());
        assert!(record.pos.is_none());
        assert!(record.mapq.is_none());
        assert_eq!(record.flag, sam::record::Flags::UNMAPPED);
        assert!(record.next_ref_id.is_none());
        assert!(record.next_pos.is_none());
        assert_eq!(record.tlen, 0);
        assert!(record.read_name.is_none());
        assert!(record.cigar.is_empty());
        assert!(record.seq.is_empty());
        assert!(record.qual.is_empty());
        assert!(record.data.is_empty());
    }
}
