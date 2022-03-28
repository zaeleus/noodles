//! BAM record and fields.

pub mod builder;
mod convert;
pub mod reference_sequence_id;

pub use self::builder::Builder;

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
    reference_sequence_id: Option<usize>,
    position: Option<Position>,
    mapping_quality: Option<sam::record::MappingQuality>,
    flags: sam::record::Flags,
    mate_reference_sequence_id: Option<usize>,
    mate_position: Option<Position>,
    template_length: i32,
    read_name: Option<sam::record::ReadName>,
    cigar: sam::record::Cigar,
    sequence: sam::record::Sequence,
    quality_scores: sam::record::QualityScores,
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
    /// let record = builder.build();
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
    pub fn reference_sequence_id(&self) -> Option<usize> {
        self.reference_sequence_id
    }

    /// Returns a mutable reference to the reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let mut record = bam::Record::default();
    /// *record.reference_sequence_id_mut() = Some(1);
    /// assert_eq!(record.reference_sequence_id(), Some(1));
    /// ```
    pub fn reference_sequence_id_mut(&mut self) -> &mut Option<usize> {
        &mut self.reference_sequence_id
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
        self.position
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
        &mut self.position
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
    /// *record.mapping_quality_mut() = MappingQuality::new(13);
    ///
    /// assert_eq!(record.mapping_quality(), MappingQuality::new(13));
    /// ```
    pub fn mapping_quality_mut(&mut self) -> &mut Option<sam::record::MappingQuality> {
        &mut self.mapping_quality
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
        &mut self.flags
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
    pub fn mate_reference_sequence_id(&self) -> Option<usize> {
        self.mate_reference_sequence_id
    }

    /// Returns a mutable reference to the mate reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let mut record = bam::Record::default();
    /// *record.mate_reference_sequence_id_mut() = Some(1);
    /// assert_eq!(record.mate_reference_sequence_id(), Some(1));
    /// ```
    pub fn mate_reference_sequence_id_mut(&mut self) -> &mut Option<usize> {
        &mut self.mate_reference_sequence_id
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
        self.mate_position
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
        &mut self.mate_position
    }

    pub(crate) fn template_length_mut(&mut self) -> &mut i32 {
        &mut self.template_length
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

    /// Returns a mutable reference to the CIGAR.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{record::Cigar, AlignmentRecord};
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
        &mut self.sequence
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
        &mut self.quality_scores
    }

    /// Returns a mutable reference to the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_sam::{
    ///     record::data::{field::{Tag, Value}, Field},
    ///     AlignmentRecord
    /// };
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
        self.flags
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
        self.mapping_quality
    }

    fn cigar(&self) -> &sam::record::Cigar {
        &self.cigar
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
        self.template_length
    }

    fn sequence(&self) -> &sam::record::Sequence {
        &self.sequence
    }

    fn quality_scores(&self) -> &sam::record::QualityScores {
        &self.quality_scores
    }

    fn data(&self) -> &sam::record::Data {
        &self.data
    }
}

fn get_reference_sequence(
    reference_sequences: &ReferenceSequences,
    reference_sequence_id: Option<usize>,
) -> Option<io::Result<&ReferenceSequence>> {
    reference_sequence_id.map(|id| {
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
        Self::builder().build()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let record = Record::default();

        assert!(record.reference_sequence_id.is_none());
        assert!(record.position.is_none());
        assert!(record.mapping_quality.is_none());
        assert_eq!(record.flags, sam::record::Flags::UNMAPPED);
        assert!(record.mate_reference_sequence_id.is_none());
        assert!(record.mate_position.is_none());
        assert_eq!(record.template_length, 0);
        assert!(record.read_name.is_none());
        assert!(record.cigar.is_empty());
        assert!(record.sequence.is_empty());
        assert!(record.quality_scores.is_empty());
        assert!(record.data.is_empty());
    }
}
