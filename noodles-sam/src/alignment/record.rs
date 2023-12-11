//! Alignment record.

mod builder;
mod flags;
mod read_name;

pub use self::builder::Builder;

use std::io;

use noodles_core::Position;

use crate::{
    header::{
        record::value::{
            map::{self, ReferenceSequence},
            Map,
        },
        ReferenceSequences,
    },
    record::{self, Cigar, Data, MappingQuality, QualityScores, Sequence},
    Header,
};

#[doc(hidden)]
pub use self::{flags::Flags, read_name::ReadName};

/// An alignment record.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    read_name: Option<record::ReadName>,
    flags: record::Flags,
    reference_sequence_id: Option<usize>,
    alignment_start: Option<Position>,
    mapping_quality: Option<MappingQuality>,
    cigar: Cigar,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<Position>,
    template_length: i32,
    sequence: Sequence,
    quality_scores: QualityScores,
    data: Data,
}

impl Record {
    /// Creates an alignment record builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let builder = sam::alignment::Record::builder();
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
    /// let record = sam::alignment::Record::default();
    /// assert!(record.read_name().is_none());
    /// ```
    pub fn read_name(&self) -> Option<&record::ReadName> {
        self.read_name.as_ref()
    }

    /// Returns a mutable reference to the read name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::ReadName};
    ///
    /// let read_name: ReadName = "r1".parse()?;
    ///
    /// let mut record = sam::alignment::Record::default();
    /// *record.read_name_mut() = Some(read_name.clone());
    ///
    /// assert_eq!(record.read_name(), Some(&read_name));
    /// Ok::<_, sam::record::read_name::ParseError>(())
    /// ```
    pub fn read_name_mut(&mut self) -> &mut Option<record::ReadName> {
        &mut self.read_name
    }

    /// Returns the flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Flags};
    /// let record = sam::alignment::Record::default();
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
    /// let mut record = sam::alignment::Record::default();
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
    /// let record = sam::alignment::Record::default();
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
    /// let mut record = sam::alignment::Record::default();
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
    /// let record = sam::alignment::Record::default();
    /// assert!(record.alignment_start().is_none());
    /// ```
    pub fn alignment_start(&self) -> Option<Position> {
        self.alignment_start
    }

    /// Returns a mutable reference to the alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam as sam;
    /// let mut record = sam::alignment::Record::default();
    /// *record.alignment_start_mut() = Some(Position::MIN);
    /// assert_eq!(record.alignment_start(), Some(Position::MIN));
    /// ```
    pub fn alignment_start_mut(&mut self) -> &mut Option<Position> {
        &mut self.alignment_start
    }

    /// Returns the mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::Record::default();
    /// assert!(record.mapping_quality().is_none());
    /// ```
    pub fn mapping_quality(&self) -> Option<MappingQuality> {
        self.mapping_quality
    }

    /// Returns a mutable reference to the mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::MappingQuality};
    /// let mut record = sam::alignment::Record::default();
    /// *record.mapping_quality_mut() = Some(MappingQuality::MIN);
    /// assert_eq!(record.mapping_quality(), Some(MappingQuality::MIN));
    /// ```
    pub fn mapping_quality_mut(&mut self) -> &mut Option<MappingQuality> {
        &mut self.mapping_quality
    }

    /// Returns the CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::Record::default();
    /// assert!(record.cigar().is_empty());
    /// ```
    pub fn cigar(&self) -> &Cigar {
        &self.cigar
    }

    /// Returns a mutable reference to the CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Cigar};
    ///
    /// let cigar: Cigar = "4M".parse()?;
    ///
    /// let mut record = sam::alignment::Record::default();
    /// *record.cigar_mut() = cigar.clone();
    ///
    /// assert_eq!(record.cigar(), &cigar);
    /// Ok::<_, sam::record::cigar::ParseError>(())
    /// ```
    pub fn cigar_mut(&mut self) -> &mut Cigar {
        &mut self.cigar
    }

    /// Returns the mate reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::Record::default();
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
    /// let mut record = sam::alignment::Record::default();
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
    /// let record = sam::alignment::Record::default();
    /// assert!(record.mate_alignment_start().is_none());
    /// ```
    pub fn mate_alignment_start(&self) -> Option<Position> {
        self.mate_alignment_start
    }

    /// Returns a mutable reference to the mate alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam as sam;
    /// let mut record = sam::alignment::Record::default();
    /// *record.mate_alignment_start_mut() = Some(Position::MIN);
    /// assert_eq!(record.mate_alignment_start(), Some(Position::MIN));
    /// ```
    pub fn mate_alignment_start_mut(&mut self) -> &mut Option<Position> {
        &mut self.mate_alignment_start
    }

    /// Returns the template length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::alignment::Record::default();
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
    /// let mut record = sam::alignment::Record::default();
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
    /// let record = sam::alignment::Record::default();
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
    /// use noodles_sam::{self as sam, record::Sequence};
    ///
    /// let sequence: Sequence = "ACGT".parse()?;
    ///
    /// let mut record = sam::alignment::Record::default();
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
    /// let record = sam::alignment::Record::default();
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
    /// use noodles_sam::{self as sam, record::QualityScores};
    ///
    /// let quality_scores: QualityScores = "NDLS".parse()?;
    ///
    /// let mut record = sam::alignment::Record::default();
    /// *record.quality_scores_mut() = quality_scores.clone();
    ///
    /// assert_eq!(record.quality_scores(), &quality_scores);
    /// Ok::<_, sam::record::quality_scores::ParseError>(())
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
    /// let record = sam::alignment::Record::default();
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
    /// use noodles_sam::{self as sam, record::Data};
    ///
    /// let data: Data = "NH:i:1".parse()?;
    ///
    /// let mut record = sam::alignment::Record::default();
    /// *record.data_mut() = data.clone();
    ///
    /// assert_eq!(record.data_mut(), &data);
    /// Ok::<_, sam::record::data::ParseError>(())
    /// ```
    pub fn data_mut(&mut self) -> &mut Data {
        &mut self.data
    }

    /// Returns the associated reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let header = sam::Header::default();
    /// let record = sam::alignment::Record::default();
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
    /// let record = sam::alignment::Record::default();
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
    /// let record = sam::alignment::Record::default();
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
    /// use noodles_sam as sam;
    ///
    /// let record = sam::alignment::Record::builder()
    ///     .set_alignment_start(Position::try_from(8)?)
    ///     .set_cigar("5M".parse()?)
    ///     .build();
    ///
    /// assert_eq!(record.alignment_end(), Position::new(12));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn alignment_end(&self) -> Option<Position> {
        self.alignment_start().and_then(|alignment_start| {
            let end = usize::from(alignment_start) + self.alignment_span() - 1;
            Position::new(end)
        })
    }
}

impl Default for Record {
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
