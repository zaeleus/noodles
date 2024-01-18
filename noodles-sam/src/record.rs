//! SAM record.

mod bounds;
mod cigar;
pub mod data;
mod flags;
mod mapping_quality;
mod position;
pub mod quality_scores;
mod read_name;
mod reference_sequence_id;
mod reference_sequence_name;
mod sequence;
mod template_length;

use std::{fmt, io};

use noodles_core as core;

use self::bounds::Bounds;
pub use self::{
    cigar::Cigar, data::Data, flags::Flags, mapping_quality::MappingQuality, position::Position,
    quality_scores::QualityScores, read_name::ReadName, reference_sequence_id::ReferenceSequenceId,
    reference_sequence_name::ReferenceSequenceName, sequence::Sequence,
    template_length::TemplateLength,
};
use crate::Header;

const MISSING: &[u8] = b"*";

/// A SAM record.
#[derive(Clone, Eq, PartialEq)]
pub struct Record {
    pub(crate) buf: Vec<u8>,
    pub(crate) bounds: Bounds,
}

impl Record {
    /// Returns the read name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::Record::default();
    /// assert!(record.read_name().is_none());
    /// ```
    pub fn read_name(&self) -> Option<ReadName<'_>> {
        match &self.buf[self.bounds.read_name_range()] {
            MISSING => None,
            buf => Some(ReadName::new(buf)),
        }
    }

    /// Returns the flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, alignment::record::Flags};
    /// let record = sam::Record::default();
    /// assert_eq!(Flags::try_from(record.flags())?, Flags::UNMAPPED);
    /// # Ok::<_, lexical_core::Error>(())
    /// ```
    pub fn flags(&self) -> Flags<'_> {
        let src = &self.buf[self.bounds.flags_range()];
        Flags::new(src)
    }

    /// Returns the reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let header = sam::Header::default();
    /// let record = sam::Record::default();
    /// assert!(record.reference_sequence_id(&header).is_none());
    /// ```
    pub fn reference_sequence_id<'r, 'h>(
        &'r self,
        header: &'h Header,
    ) -> Option<ReferenceSequenceId<'h, 'r>> {
        self.reference_sequence_name()
            .map(|reference_sequence_name| {
                ReferenceSequenceId::new(header, reference_sequence_name)
            })
    }

    /// Returns the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::Record::default();
    /// assert!(record.reference_sequence_name().is_none());
    /// ```
    pub fn reference_sequence_name(&self) -> Option<ReferenceSequenceName<'_>> {
        match &self.buf[self.bounds.reference_sequence_name_range()] {
            MISSING => None,
            buf => Some(ReferenceSequenceName::new(buf)),
        }
    }

    /// Returns the alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::Record::default();
    /// assert!(record.alignment_start().is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn alignment_start(&self) -> Option<Position<'_>> {
        const MISSING: &[u8] = b"0";

        match &self.buf[self.bounds.alignment_start_range()] {
            MISSING => None,
            buf => Some(Position::new(buf)),
        }
    }

    /// Returns the mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::Record::default();
    /// assert!(record.mapping_quality().is_none());
    /// ```
    pub fn mapping_quality(&self) -> Option<MappingQuality<'_>> {
        const MISSING: &[u8] = b"255";

        match &self.buf[self.bounds.mapping_quality_range()] {
            MISSING => None,
            buf => Some(MappingQuality::new(buf)),
        }
    }

    /// Returns the CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::Record::default();
    /// assert!(record.cigar().is_empty());
    /// ```
    pub fn cigar(&self) -> Cigar<'_> {
        match &self.buf[self.bounds.cigar_range()] {
            MISSING => Cigar::new(b""),
            buf => Cigar::new(buf),
        }
    }

    /// Returns the mate reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let header = sam::Header::default();
    /// let record = sam::Record::default();
    /// assert!(record.mate_reference_sequence_id(&header).is_none());
    /// ```
    pub fn mate_reference_sequence_id<'r, 'h>(
        &'r self,
        header: &'h Header,
    ) -> Option<ReferenceSequenceId<'h, 'r>> {
        self.mate_reference_sequence_name()
            .map(|mate_reference_sequence_name| {
                ReferenceSequenceId::new(header, mate_reference_sequence_name)
            })
    }

    /// Returns the mate reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::Record::default();
    /// assert!(record.mate_reference_sequence_name().is_none());
    /// ```
    pub fn mate_reference_sequence_name(&self) -> Option<ReferenceSequenceName<'_>> {
        const EQ: &[u8] = b"=";

        match &self.buf[self.bounds.mate_reference_sequence_name_range()] {
            MISSING => None,
            EQ => self.reference_sequence_name(),
            buf => Some(ReferenceSequenceName::new(buf)),
        }
    }

    /// Returns the mate alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::Record::default();
    /// assert!(record.mate_alignment_start().is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn mate_alignment_start(&self) -> Option<Position<'_>> {
        const MISSING: &[u8] = b"0";

        match &self.buf[self.bounds.mate_alignment_start_range()] {
            MISSING => None,
            buf => Some(Position::new(buf)),
        }
    }

    /// Returns the template length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::Record::default();
    /// assert_eq!(i32::try_from(record.template_length())?, 0);
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn template_length(&self) -> TemplateLength<'_> {
        let buf = &self.buf[self.bounds.template_length_range()];
        TemplateLength::new(buf)
    }

    /// Returns the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::Record::default();
    /// assert!(record.sequence().is_empty());
    /// ```
    pub fn sequence(&self) -> Sequence<'_> {
        let buf = match &self.buf[self.bounds.sequence_range()] {
            MISSING => b"",
            buf => buf,
        };

        Sequence::new(buf)
    }

    /// Returns the quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::Record::default();
    /// assert!(record.quality_scores().is_empty());
    /// ```
    pub fn quality_scores(&self) -> QualityScores<'_> {
        let buf = match &self.buf[self.bounds.quality_scores_range()] {
            MISSING => b"",
            buf => buf,
        };

        QualityScores::new(buf)
    }

    /// Returns the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::Record::default();
    /// assert!(record.data().is_empty());
    /// ```
    pub fn data(&self) -> Data<'_> {
        let buf = &self.buf[self.bounds.data_range()];
        Data::new(buf)
    }
}

impl fmt::Debug for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Record")
            .field("read_name", &self.read_name())
            .field("flags", &self.flags())
            .field("reference_sequence_name", &self.reference_sequence_name())
            .field("alignment_start", &self.alignment_start())
            .field("mapping_quality", &self.mapping_quality())
            .field("cigar", &self.cigar())
            .field(
                "mate_reference_sequence_name",
                &self.mate_reference_sequence_name(),
            )
            .field("mate_alignment_start", &self.mate_alignment_start())
            .field("template_length", &self.template_length())
            .field("sequence", &self.sequence())
            .field("quality_scores", &self.quality_scores())
            .field("data", &self.data())
            .finish()
    }
}

impl crate::alignment::Record for Record {
    fn name(&self) -> Option<Box<dyn crate::alignment::record::fields::Name + '_>> {
        let read_name = self.read_name()?;
        Some(Box::new(read_name))
    }

    fn flags(&self) -> io::Result<crate::alignment::record::Flags> {
        crate::alignment::record::Flags::try_from(self.flags())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    fn reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        header: &'h Header,
    ) -> Option<io::Result<usize>> {
        self.reference_sequence_id(header).map(usize::try_from)
    }

    fn alignment_start(&self) -> Option<io::Result<core::Position>> {
        self.alignment_start().map(core::Position::try_from)
    }

    fn mapping_quality(&self) -> Option<io::Result<crate::alignment::record::MappingQuality>> {
        self.mapping_quality()
            .map(crate::alignment::record::MappingQuality::try_from)
    }

    fn cigar(&self) -> Box<dyn crate::alignment::record::fields::Cigar + '_> {
        Box::new(self.cigar())
    }

    fn mate_reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        header: &'h Header,
    ) -> Option<io::Result<usize>> {
        self.mate_reference_sequence_id(header).map(usize::try_from)
    }

    fn mate_alignment_start(&self) -> Option<io::Result<core::Position>> {
        self.mate_alignment_start().map(core::Position::try_from)
    }

    fn template_length(&self) -> io::Result<i32> {
        i32::try_from(self.template_length())
    }

    fn sequence(&self) -> Box<dyn crate::alignment::record::fields::Sequence + '_> {
        Box::new(self.sequence())
    }

    fn quality_scores(&self) -> Box<dyn crate::alignment::record::fields::QualityScores + '_> {
        Box::new(self.quality_scores())
    }

    fn data(&self) -> Box<dyn crate::alignment::record::fields::Data + '_> {
        Box::new(self.data())
    }
}

impl Default for Record {
    fn default() -> Self {
        let buf = b"*4*0255**00**".to_vec();

        let bounds = Bounds {
            read_name_end: 1,
            flags_end: 2,
            reference_sequence_name_end: 3,
            alignment_start_end: 4,
            mapping_quality_end: 7,
            cigar_end: 8,
            mate_reference_sequence_name_end: 9,
            mate_alignment_start_end: 10,
            template_length_end: 11,
            sequence_end: 12,
            quality_scores_end: 13,
        };

        Self { buf, bounds }
    }
}
