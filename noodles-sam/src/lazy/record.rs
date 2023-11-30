mod bounds;
mod cigar;
mod data;
mod flags;
mod position;
mod quality_scores;
mod reference_sequence_name;
mod sequence;

use std::{fmt, io};

use self::bounds::Bounds;
pub use self::{
    cigar::Cigar, data::Data, flags::Flags, position::Position, quality_scores::QualityScores,
    reference_sequence_name::ReferenceSequenceName, sequence::Sequence,
};
use crate::record::MappingQuality;

const MISSING: &[u8] = b"*";

/// An immutable, lazily-evalulated SAM record.
///
/// The fields are _not_ memoized.
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
    /// let record = sam::lazy::Record::default();
    /// assert!(record.read_name().is_none());
    /// ```
    pub fn read_name(&self) -> Option<&[u8]> {
        match &self.buf[self.bounds.read_name_range()] {
            MISSING => None,
            buf => Some(buf),
        }
    }

    /// Returns the flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Flags};
    /// let record = sam::lazy::Record::default();
    /// assert_eq!(Flags::try_from(record.flags())?, Flags::UNMAPPED);
    /// # Ok::<_, lexical_core::Error>(())
    /// ```
    pub fn flags(&self) -> Flags<'_> {
        let src = &self.buf[self.bounds.flags_range()];
        Flags::new(src)
    }

    /// Returns the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::lazy::Record::default();
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
    /// let record = sam::lazy::Record::default();
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
    /// let record = sam::lazy::Record::default();
    /// assert!(record.mapping_quality()?.is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn mapping_quality(&self) -> io::Result<Option<MappingQuality>> {
        use crate::reader::record::parse_mapping_quality;
        let src = &self.buf[self.bounds.mapping_quality_range()];
        parse_mapping_quality(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    /// Returns the CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::lazy::Record::default();
    /// assert!(record.cigar().is_empty());
    /// ```
    pub fn cigar(&self) -> Cigar<'_> {
        match &self.buf[self.bounds.cigar_range()] {
            MISSING => Cigar::new(b""),
            buf => Cigar::new(buf),
        }
    }

    /// Returns the mate reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::lazy::Record::default();
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
    /// let record = sam::lazy::Record::default();
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
    /// let record = sam::lazy::Record::default();
    /// assert_eq!(record.template_length()?, 0);
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn template_length(&self) -> io::Result<i32> {
        use crate::reader::record::parse_template_length;
        let src = &self.buf[self.bounds.template_length_range()];
        parse_template_length(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    /// Returns the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::lazy::Record::default();
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
    /// let record = sam::lazy::Record::default();
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
    /// let record = sam::lazy::Record::default();
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
