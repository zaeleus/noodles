//! Alignment record buffer.

mod builder;
mod name;
mod quality_scores;
mod sequence;

pub use self::{builder::Builder, name::Name, quality_scores::QualityScores, sequence::Sequence};

use std::io;

use noodles_core as core;

use super::Record;
use crate::{
    header::{
        record::value::{
            map::{self, ReferenceSequence},
            Map,
        },
        ReferenceSequences,
    },
    record, Header,
};

/// An alignment record.
#[derive(Clone, Debug, PartialEq)]
pub struct RecordBuf {
    name: Option<Name>,
    flags: record::Flags,
    reference_sequence_id: Option<usize>,
    alignment_start: Option<core::Position>,
    mapping_quality: Option<record::MappingQuality>,
    cigar: record::Cigar,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<core::Position>,
    template_length: i32,
    sequence: Sequence,
    quality_scores: QualityScores,
    data: record::Data,
}

impl RecordBuf {
    /// Creates an alignment record builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let builder = sam::alignment::RecordBuf::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns the read name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::RecordBuf::default();
    /// assert!(record.name().is_none());
    /// ```
    pub fn name(&self) -> Option<&Name> {
        self.name.as_ref()
    }

    /// Returns a mutable reference to the read name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, alignment::record_buf::Name};
    ///
    /// let name = Name::from(b"r1");
    ///
    /// let mut record = sam::alignment::RecordBuf::default();
    /// *record.name_mut() = Some(name.clone());
    ///
    /// assert_eq!(record.name(), Some(&name));
    /// ```
    pub fn name_mut(&mut self) -> &mut Option<Name> {
        &mut self.name
    }

    /// Returns the flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Flags};
    /// let record = sam::alignment::RecordBuf::default();
    /// assert_eq!(record.flags(), Flags::UNMAPPED);
    /// ```
    pub fn flags(&self) -> record::Flags {
        self.flags
    }

    /// Returns a mutable reference to the flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Flags};
    /// let mut record = sam::alignment::RecordBuf::default();
    /// *record.flags_mut() = Flags::empty();
    /// assert!(record.flags().is_empty());
    /// ```
    pub fn flags_mut(&mut self) -> &mut record::Flags {
        &mut self.flags
    }

    /// Returns the reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::RecordBuf::default();
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
    /// use noodles_sam as sam;
    /// let mut record = sam::alignment::RecordBuf::default();
    /// *record.reference_sequence_id_mut() = Some(0);
    /// assert_eq!(record.reference_sequence_id(), Some(0));
    /// ```
    pub fn reference_sequence_id_mut(&mut self) -> &mut Option<usize> {
        &mut self.reference_sequence_id
    }

    /// Returns the alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::RecordBuf::default();
    /// assert!(record.alignment_start().is_none());
    /// ```
    pub fn alignment_start(&self) -> Option<core::Position> {
        self.alignment_start
    }

    /// Returns a mutable reference to the alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam as sam;
    /// let mut record = sam::alignment::RecordBuf::default();
    /// *record.alignment_start_mut() = Some(Position::MIN);
    /// assert_eq!(record.alignment_start(), Some(Position::MIN));
    /// ```
    pub fn alignment_start_mut(&mut self) -> &mut Option<core::Position> {
        &mut self.alignment_start
    }

    /// Returns the mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::RecordBuf::default();
    /// assert!(record.mapping_quality().is_none());
    /// ```
    pub fn mapping_quality(&self) -> Option<record::MappingQuality> {
        self.mapping_quality
    }

    /// Returns a mutable reference to the mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::MappingQuality};
    /// let mut record = sam::alignment::RecordBuf::default();
    /// *record.mapping_quality_mut() = Some(MappingQuality::MIN);
    /// assert_eq!(record.mapping_quality(), Some(MappingQuality::MIN));
    /// ```
    pub fn mapping_quality_mut(&mut self) -> &mut Option<record::MappingQuality> {
        &mut self.mapping_quality
    }

    /// Returns the CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::RecordBuf::default();
    /// assert!(record.cigar().is_empty());
    /// ```
    pub fn cigar(&self) -> &record::Cigar {
        &self.cigar
    }

    /// Returns a mutable reference to the CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     self as sam,
    ///     record::{
    ///         cigar::{op::Kind, Op},
    ///         Cigar,
    ///     },
    /// };
    ///
    /// let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
    ///
    /// let mut record = sam::alignment::RecordBuf::default();
    /// *record.cigar_mut() = cigar.clone();
    ///
    /// assert_eq!(record.cigar(), &cigar);
    /// ```
    pub fn cigar_mut(&mut self) -> &mut record::Cigar {
        &mut self.cigar
    }

    /// Returns the mate reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::RecordBuf::default();
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
    /// use noodles_sam as sam;
    /// let mut record = sam::alignment::RecordBuf::default();
    /// *record.mate_reference_sequence_id_mut() = Some(0);
    /// assert_eq!(record.mate_reference_sequence_id(), Some(0));
    /// ```
    pub fn mate_reference_sequence_id_mut(&mut self) -> &mut Option<usize> {
        &mut self.mate_reference_sequence_id
    }

    /// Returns the mate alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::RecordBuf::default();
    /// assert!(record.mate_alignment_start().is_none());
    /// ```
    pub fn mate_alignment_start(&self) -> Option<core::Position> {
        self.mate_alignment_start
    }

    /// Returns a mutable reference to the mate alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam as sam;
    /// let mut record = sam::alignment::RecordBuf::default();
    /// *record.mate_alignment_start_mut() = Some(Position::MIN);
    /// assert_eq!(record.mate_alignment_start(), Some(Position::MIN));
    /// ```
    pub fn mate_alignment_start_mut(&mut self) -> &mut Option<core::Position> {
        &mut self.mate_alignment_start
    }

    /// Returns the template length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::RecordBuf::default();
    /// assert_eq!(record.template_length(), 0);
    /// ```
    pub fn template_length(&self) -> i32 {
        self.template_length
    }

    /// Returns a mutable reference to template length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let mut record = sam::alignment::RecordBuf::default();
    /// *record.template_length_mut() = 4;
    /// assert_eq!(record.template_length(), 4);
    /// ```
    pub fn template_length_mut(&mut self) -> &mut i32 {
        &mut self.template_length
    }

    /// Returns the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::RecordBuf::default();
    /// assert!(record.sequence().is_empty());
    /// ```
    pub fn sequence(&self) -> &Sequence {
        &self.sequence
    }

    /// Returns a mutable reference to sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, alignment::record_buf::Sequence};
    ///
    /// let sequence = Sequence::from(b"ACGT".to_vec());
    ///
    /// let mut record = sam::alignment::RecordBuf::default();
    /// *record.sequence_mut() = sequence.clone();
    ///
    /// assert_eq!(record.sequence(), &sequence);
    /// Ok::<_, sam::record::sequence::ParseError>(())
    /// ```
    pub fn sequence_mut(&mut self) -> &mut Sequence {
        &mut self.sequence
    }

    /// Returns the quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::RecordBuf::default();
    /// assert!(record.quality_scores().is_empty());
    /// ```
    pub fn quality_scores(&self) -> &QualityScores {
        &self.quality_scores
    }

    /// Returns a mutable reference to quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, alignment::record_buf::QualityScores};
    ///
    /// let quality_scores = QualityScores::default();
    ///
    /// let mut record = sam::alignment::RecordBuf::default();
    /// *record.quality_scores_mut() = quality_scores.clone();
    ///
    /// assert_eq!(record.quality_scores(), &quality_scores);
    /// ```
    pub fn quality_scores_mut(&mut self) -> &mut QualityScores {
        &mut self.quality_scores
    }

    /// Returns the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::RecordBuf::default();
    /// assert!(record.data().is_empty());
    /// ```
    pub fn data(&self) -> &record::Data {
        &self.data
    }

    /// Returns a mutable reference to the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Data};
    ///
    /// let data: Data = "NH:i:1".parse()?;
    ///
    /// let mut record = sam::alignment::RecordBuf::default();
    /// *record.data_mut() = data.clone();
    ///
    /// assert_eq!(record.data_mut(), &data);
    /// Ok::<_, sam::record::data::ParseError>(())
    /// ```
    pub fn data_mut(&mut self) -> &mut record::Data {
        &mut self.data
    }

    /// Returns the associated reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let header = sam::Header::default();
    /// let record = sam::alignment::RecordBuf::default();
    /// assert!(record.reference_sequence(&header).is_none());
    /// ```
    pub fn reference_sequence<'a>(
        &self,
        header: &'a Header,
    ) -> Option<
        io::Result<(
            &'a map::reference_sequence::Name,
            &'a Map<ReferenceSequence>,
        )>,
    > {
        get_reference_sequence(header.reference_sequences(), self.reference_sequence_id())
    }

    /// Returns the associated mate reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let header = sam::Header::default();
    /// let record = sam::alignment::RecordBuf::default();
    /// assert!(record.mate_reference_sequence(&header).is_none());
    /// ```
    pub fn mate_reference_sequence<'a>(
        &self,
        header: &'a Header,
    ) -> Option<
        io::Result<(
            &'a map::reference_sequence::Name,
            &'a Map<ReferenceSequence>,
        )>,
    > {
        get_reference_sequence(
            header.reference_sequences(),
            self.mate_reference_sequence_id(),
        )
    }

    /// Returns the alignment span.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::RecordBuf::default();
    /// assert_eq!(record.alignment_span(), 0);
    /// ```
    pub fn alignment_span(&self) -> usize {
        self.cigar().alignment_span()
    }

    /// Calculates the end position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam::{
    ///     self as sam,
    ///     record::cigar::{op::Kind, Op},
    /// };
    ///
    /// let record = sam::alignment::RecordBuf::builder()
    ///     .set_alignment_start(Position::try_from(8)?)
    ///     .set_cigar([Op::new(Kind::Match, 5)].into_iter().collect())
    ///     .build();
    ///
    /// assert_eq!(record.alignment_end(), Position::new(12));
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn alignment_end(&self) -> Option<core::Position> {
        self.alignment_start().and_then(|alignment_start| {
            let end = usize::from(alignment_start) + self.alignment_span() - 1;
            core::Position::new(end)
        })
    }
}

impl Record for RecordBuf {
    fn name(&self) -> Option<Box<dyn super::record::Name + '_>> {
        let name = self.name()?;
        Some(Box::new(name))
    }

    fn flags(&self) -> Box<dyn super::record::Flags + '_> {
        Box::new(self.flags())
    }

    fn reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        _: &'h Header,
    ) -> Option<Box<dyn super::record::ReferenceSequenceId + 'r>> {
        let reference_sequence_id = self.reference_sequence_id()?;
        Some(Box::new(reference_sequence_id))
    }

    fn alignment_start(&self) -> Option<Box<dyn super::record::Position + '_>> {
        let alignment_start = self.alignment_start()?;
        Some(Box::new(alignment_start))
    }

    fn mapping_quality(&self) -> Option<Box<dyn super::record::MappingQuality + '_>> {
        let mapping_quality = self.mapping_quality()?;
        Some(Box::new(mapping_quality))
    }

    fn cigar(&self, _: &Header) -> Box<dyn super::record::Cigar + '_> {
        Box::new(self.cigar())
    }

    fn mate_reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        _: &'h Header,
    ) -> Option<Box<dyn super::record::ReferenceSequenceId + 'r>> {
        let mate_reference_sequence_id = self.mate_reference_sequence_id()?;
        Some(Box::new(mate_reference_sequence_id))
    }

    fn mate_alignment_start(&self) -> Option<Box<dyn super::record::Position + '_>> {
        let mate_alignment_start = self.mate_alignment_start()?;
        Some(Box::new(mate_alignment_start))
    }

    fn template_length(&self) -> Box<dyn super::record::TemplateLength + '_> {
        Box::new(self.template_length())
    }

    fn sequence(&self) -> Box<dyn super::record::Sequence + '_> {
        Box::new(self.sequence())
    }

    fn quality_scores(&self) -> Box<dyn super::record::QualityScores + '_> {
        Box::new(self.quality_scores())
    }

    fn data(&self) -> Box<dyn super::record::Data + '_> {
        Box::new(self.data())
    }
}

impl Default for RecordBuf {
    fn default() -> Self {
        Self::builder().build()
    }
}

fn get_reference_sequence(
    reference_sequences: &ReferenceSequences,
    reference_sequence_id: Option<usize>,
) -> Option<io::Result<(&map::reference_sequence::Name, &Map<ReferenceSequence>)>> {
    reference_sequence_id.map(|id| {
        reference_sequences.get_index(id).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "invalid reference sequence ID")
        })
    })
}
