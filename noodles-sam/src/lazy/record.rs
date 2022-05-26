use std::{
    fmt, io,
    ops::{Range, RangeFrom},
    str,
};

use noodles_core::Position;

use crate::record::{
    Cigar, Data, Flags, MappingQuality, QualityScores, ReadName, ReferenceSequenceName, Sequence,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Bounds {
    pub(crate) read_name_end: usize,
    pub(crate) flags_end: usize,
    pub(crate) reference_sequence_name_end: usize,
    pub(crate) alignment_start_end: usize,
    pub(crate) mapping_quality_end: usize,
    pub(crate) cigar_end: usize,
    pub(crate) mate_reference_sequence_name_end: usize,
    pub(crate) mate_alignment_start_end: usize,
    pub(crate) template_length_end: usize,
    pub(crate) sequence_end: usize,
    pub(crate) quality_scores_end: usize,
}

impl Bounds {
    fn read_name_range(&self) -> Range<usize> {
        0..self.read_name_end
    }

    fn flags_range(&self) -> Range<usize> {
        self.read_name_end..self.flags_end
    }

    fn reference_sequence_name_range(&self) -> Range<usize> {
        self.flags_end..self.reference_sequence_name_end
    }

    fn alignment_start_range(&self) -> Range<usize> {
        self.reference_sequence_name_end..self.alignment_start_end
    }

    fn mapping_quality_range(&self) -> Range<usize> {
        self.alignment_start_end..self.mapping_quality_end
    }

    fn cigar_range(&self) -> Range<usize> {
        self.mapping_quality_end..self.cigar_end
    }

    fn mate_reference_sequence_name_range(&self) -> Range<usize> {
        self.cigar_end..self.mate_reference_sequence_name_end
    }

    fn mate_alignment_start_range(&self) -> Range<usize> {
        self.mate_reference_sequence_name_end..self.mate_alignment_start_end
    }

    fn template_length_range(&self) -> Range<usize> {
        self.mate_alignment_start_end..self.template_length_end
    }

    fn sequence_range(&self) -> Range<usize> {
        self.template_length_end..self.sequence_end
    }

    fn quality_scores_range(&self) -> Range<usize> {
        self.sequence_end..self.quality_scores_end
    }

    fn data_range(&self) -> RangeFrom<usize> {
        self.quality_scores_end..
    }
}

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
    /// assert!(record.read_name()?.is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_name(&self) -> io::Result<Option<ReadName>> {
        const MISSING: &[u8] = b"*";

        let buf = &self.buf[self.bounds.read_name_range()];

        match buf {
            MISSING => Ok(None),
            _ => ReadName::try_new(buf)
                .map(Some)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        }
    }

    /// Returns the flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Flags};
    /// let record = sam::lazy::Record::default();
    /// assert_eq!(record.flags()?, Flags::UNMAPPED);
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn flags(&self) -> io::Result<Flags> {
        let buf = &self.buf[self.bounds.flags_range()];

        str::from_utf8(buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                s.parse::<u16>()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(Flags::from)
    }

    /// Returns the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::lazy::Record::default();
    /// assert!(record.reference_sequence_name()?.is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn reference_sequence_name(&self) -> io::Result<Option<ReferenceSequenceName>> {
        let buf = &self.buf[self.bounds.reference_sequence_name_range()];
        parse_reference_sequence_name(buf)
    }

    /// Returns the alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::lazy::Record::default();
    /// assert!(record.alignment_start()?.is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn alignment_start(&self) -> io::Result<Option<Position>> {
        let buf = &self.buf[self.bounds.alignment_start_range()];
        parse_alignment_start(buf)
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
        let buf = &self.buf[self.bounds.mapping_quality_range()];

        str::from_utf8(buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                s.parse()
                    .map(MappingQuality::new)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    /// Returns the CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::lazy::Record::default();
    /// assert!(record.cigar()?.is_empty());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn cigar(&self) -> io::Result<Cigar> {
        let buf = &self.buf[self.bounds.cigar_range()];

        str::from_utf8(buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| match s {
                "*" => Ok(Cigar::default()),
                _ => s
                    .parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            })
    }

    /// Returns the mate reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::lazy::Record::default();
    /// assert!(record.mate_reference_sequence_name()?.is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn mate_reference_sequence_name(&self) -> io::Result<Option<ReferenceSequenceName>> {
        let buf = &self.buf[self.bounds.mate_reference_sequence_name_range()];
        parse_reference_sequence_name(buf)
    }

    /// Returns the mate alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::lazy::Record::default();
    /// assert!(record.mate_alignment_start()?.is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn mate_alignment_start(&self) -> io::Result<Option<Position>> {
        let buf = &self.buf[self.bounds.mate_alignment_start_range()];
        parse_alignment_start(buf)
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
        let buf = &self.buf[self.bounds.template_length_range()];

        str::from_utf8(buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    /// Returns the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::lazy::Record::default();
    /// assert!(record.sequence()?.is_empty());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn sequence(&self) -> io::Result<Sequence> {
        let buf = &self.buf[self.bounds.sequence_range()];

        str::from_utf8(buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| match s {
                "*" => Ok(Sequence::default()),
                _ => s
                    .parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            })
    }

    /// Returns the quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::lazy::Record::default();
    /// assert!(record.quality_scores()?.is_empty());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn quality_scores(&self) -> io::Result<QualityScores> {
        let buf = &self.buf[self.bounds.quality_scores_range()];

        str::from_utf8(buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| match s {
                "*" => Ok(QualityScores::default()),
                _ => s
                    .parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            })
    }

    /// Returns the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let record = sam::lazy::Record::default();
    /// assert!(record.data()?.is_empty());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn data(&self) -> io::Result<Data> {
        let buf = &self.buf[self.bounds.data_range()];

        str::from_utf8(buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
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

fn parse_reference_sequence_name(buf: &[u8]) -> io::Result<Option<ReferenceSequenceName>> {
    str::from_utf8(buf)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|s| match s {
            "*" => Ok(None),
            _ => s
                .parse()
                .map(Some)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        })
}

fn parse_alignment_start(buf: &[u8]) -> io::Result<Option<Position>> {
    str::from_utf8(buf)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|s| {
            s.parse()
                .map(Position::new)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
}
