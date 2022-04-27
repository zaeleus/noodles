//! SAM record and fields.

pub mod builder;
pub mod cigar;
pub mod data;
mod field;
mod flags;
pub mod mapping_quality;
mod parser;
pub mod quality_scores;
pub mod read_name;
pub mod reference_sequence_name;

pub use self::{
    builder::Builder, cigar::Cigar, data::Data, field::Field, flags::Flags,
    mapping_quality::MappingQuality, parser::ParseError, quality_scores::QualityScores,
    read_name::ReadName, reference_sequence_name::ReferenceSequenceName,
};

#[deprecated(
    since = "0.15.0",
    note = "Use `noodles_sam::alignment::record::sequence` instead."
)]
pub use super::alignment::record::sequence;

#[deprecated(
    since = "0.15.0",
    note = "Use `noodles_sam::alignment::record::Sequence` instead."
)]
// pub use super::alignment::record::Sequence;
use super::alignment::record::Sequence;

use std::{fmt, io, str::FromStr};

use noodles_core::Position;

use super::{
    header::{ReferenceSequence, ReferenceSequences},
    AlignmentRecord,
};

pub(crate) const NULL_FIELD: &str = "*";
const EQ_FIELD: &str = "=";

/// A SAM record.
///
/// A SAM record has 11 required fields:
///
///   1. read name (`QNAME`),
///   2. flags (`FLAG`),
///   3. reference sequence name (`RNAME`),
///   4. position (`POS`),
///   5. mapping quality (`MAPQ`),
///   6. CIGAR string (`CIGAR`),
///   7. mate reference sequence name (`RNEXT`),
///   8. mate position (`PNEXT`),
///   9. template length (`TLEN`),
///   10. sequence (`SEQ`), and
///   11. quality scores (`QUAL`).
///
/// Additionally, optional data fields can be included with any record.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    read_name: Option<ReadName>,
    flags: Flags,
    reference_sequence_name: Option<ReferenceSequenceName>,
    position: Option<Position>,
    mapping_quality: Option<MappingQuality>,
    cigar: Cigar,
    mate_reference_sequence_name: Option<ReferenceSequenceName>,
    mate_position: Option<Position>,
    template_length: i32,
    sequence: Sequence,
    quality_scores: QualityScores,
    data: Data,
}

impl Record {
    /// Returns a builder to create a record from each of its fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Flags, AlignmentRecord};
    ///
    /// let record = sam::Record::builder()
    ///     .set_read_name("r0".parse()?)
    ///     .set_flags(Flags::UNMAPPED)
    ///     .build();
    ///
    /// assert_eq!(record.read_name().map(|name| name.as_ref()), Some("r0"));
    /// assert_eq!(record.flags(), Flags::UNMAPPED);
    /// assert!(record.reference_sequence_name().is_none());
    /// assert!(record.position().is_none());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns a mutable reference to the read name.
    ///
    /// This is also called the query name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, AlignmentRecord};
    ///
    /// let mut record = sam::Record::default();
    /// record.read_name_mut().insert("r0".parse()?);
    /// assert_eq!(record.read_name().map(|name| name.as_ref()), Some("r0"));
    ///
    /// *record.read_name_mut() = None;
    /// assert!(record.read_name().is_none());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn read_name_mut(&mut self) -> &mut Option<ReadName> {
        &mut self.read_name
    }

    /// Returns a mutable reference to the SAM flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Flags, AlignmentRecord};
    ///
    /// let mut record = sam::Record::builder()
    ///     .set_flags(Flags::PAIRED | Flags::READ_1)
    ///     .build();
    /// record.flags_mut().set(Flags::DUPLICATE, true);
    /// assert_eq!(record.flags(), Flags::PAIRED | Flags::READ_1 | Flags::DUPLICATE);
    ///
    /// record.flags_mut().set(Flags::PAIRED | Flags::QC_FAIL, false);
    /// assert_eq!(record.flags(), Flags::READ_1 | Flags::DUPLICATE);
    /// ```
    pub fn flags_mut(&mut self) -> &mut Flags {
        &mut self.flags
    }

    /// Returns the reference sequence name of this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let record = sam::Record::default();
    /// assert!(record.reference_sequence_name().is_none());
    ///
    /// let record = sam::Record::builder()
    ///     .set_reference_sequence_name("sq0".parse()?)
    ///     .build();
    ///
    /// assert_eq!(
    ///     record.reference_sequence_name().map(|name| name.as_str()),
    ///     Some("sq0")
    /// );
    /// # Ok::<(), sam::record::reference_sequence_name::ParseError>(())
    /// ```
    pub fn reference_sequence_name(&self) -> Option<&ReferenceSequenceName> {
        self.reference_sequence_name.as_ref()
    }

    /// Returns a mutable reference to the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let mut record = sam::Record::default();
    /// *record.reference_sequence_name_mut() = Some("sq0".parse()?);
    /// assert_eq!(record.reference_sequence_name().map(|name| name.as_str()), Some("sq0"));
    ///
    /// *record.reference_sequence_name_mut() = None;
    /// assert!(record.reference_sequence_name().is_none());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_sequence_name_mut(&mut self) -> &mut Option<ReferenceSequenceName> {
        &mut self.reference_sequence_name
    }

    /// Returns the start position of this record.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam as sam;
    ///
    /// let record = sam::Record::default();
    /// assert!(record.position().is_none());
    ///
    /// let record = sam::Record::builder()
    ///     .set_position(Position::try_from(13)?)
    ///     .build();
    /// assert_eq!(record.position(), Position::new(13));
    /// # Ok::<(), noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn position(&self) -> Option<Position> {
        self.position
    }

    /// Returns a mutable reference to the start position.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// use noodles_core::Position;
    ///
    /// let mut record = sam::Record::default();
    /// *record.position_mut() = Position::new(13);
    /// assert_eq!(record.position(), Position::new(13));
    ///
    /// *record.position_mut() = None;
    /// assert!(record.position().is_none());
    /// ```
    pub fn position_mut(&mut self) -> &mut Option<Position> {
        &mut self.position
    }

    /// Returns a mutable reference to the mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::MappingQuality, AlignmentRecord};
    ///
    /// let mut record = sam::Record::default();
    /// *record.mapping_quality_mut() = MappingQuality::new(8);
    /// assert_eq!(record.mapping_quality(), MappingQuality::new(8));
    ///
    /// *record.mapping_quality_mut() = None;
    /// assert!(record.mapping_quality().is_none());
    /// ```
    pub fn mapping_quality_mut(&mut self) -> &mut Option<MappingQuality> {
        &mut self.mapping_quality
    }

    /// Returns a mutable reference to the CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     self as sam,
    ///     record::{cigar::{op::Kind, Op}, Cigar},
    ///     AlignmentRecord,
    /// };
    ///
    /// let mut record = sam::Record::default();
    /// assert!(record.cigar().is_empty());
    ///
    /// let cigar = Cigar::from(vec![
    ///     Op::new(Kind::Match, 36),
    ///     Op::new(Kind::SoftClip, 2),
    /// ]);
    /// *record.cigar_mut() = cigar.clone();
    ///
    /// assert_eq!(record.cigar(), &cigar);
    /// ```
    pub fn cigar_mut(&mut self) -> &mut Cigar {
        &mut self.cigar
    }

    /// Returns the mate reference sequence name of this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let record = sam::Record::default();
    /// assert!(record.mate_reference_sequence_name().is_none());
    ///
    /// let record = sam::Record::builder()
    ///     .set_mate_reference_sequence_name("sq0".parse()?)
    ///     .build();
    ///
    /// assert_eq!(
    ///     record.mate_reference_sequence_name().map(|name| name.as_str()),
    ///     Some("sq0")
    /// );
    /// # Ok::<(), sam::record::reference_sequence_name::ParseError>(())
    /// ```
    pub fn mate_reference_sequence_name(&self) -> Option<&ReferenceSequenceName> {
        self.mate_reference_sequence_name.as_ref()
    }

    /// Returns a mutable reference to the mate reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let mut record = sam::Record::default();
    /// *record.mate_reference_sequence_name_mut() = Some("sq0".parse()?);
    /// assert_eq!(record.mate_reference_sequence_name().map(|name| name.as_str()), Some("sq0"));
    ///
    /// *record.mate_reference_sequence_name_mut() = None;
    /// assert!(record.mate_reference_sequence_name().is_none());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn mate_reference_sequence_name_mut(&mut self) -> &mut Option<ReferenceSequenceName> {
        &mut self.mate_reference_sequence_name
    }

    /// Returns the start position of the mate of this record.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam as sam;
    ///
    /// let record = sam::Record::default();
    /// assert!(record.mate_position().is_none());
    ///
    /// let record = sam::Record::builder()
    ///     .set_mate_position(Position::try_from(21)?)
    ///     .build();
    /// assert_eq!(record.mate_position(), Position::new(21));
    /// # Ok::<(), noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn mate_position(&self) -> Option<Position> {
        self.mate_position
    }

    /// Returns a mutable reference to the start position of the mate.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam as sam;
    ///
    /// let mut record = sam::Record::default();
    /// *record.mate_position_mut() = Position::new(13);
    /// assert_eq!(record.mate_position(), Position::new(13));
    ///
    /// *record.mate_position_mut() = None;
    /// assert!(record.mate_position().is_none());
    /// ```
    pub fn mate_position_mut(&mut self) -> &mut Option<Position> {
        &mut self.mate_position
    }

    /// Returns a mutable reference to the template length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::MappingQuality, AlignmentRecord};
    ///
    /// let mut record = sam::Record::default();
    /// *record.template_length_mut() = 101;
    /// assert_eq!(record.template_length(), 101);
    /// ```
    pub fn template_length_mut(&mut self) -> &mut i32 {
        &mut self.template_length
    }

    /// Returns a mutable reference to the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, alignment::record::Sequence, AlignmentRecord};
    ///
    /// let mut record = sam::Record::default();
    /// assert!(record.sequence().is_empty());
    ///
    /// let sequence: Sequence = "ACGT".parse()?;
    /// *record.sequence_mut() = sequence.clone();
    ///
    /// assert_eq!(record.sequence(), sequence);
    /// # Ok::<(), sam::alignment::record::sequence::ParseError>(())
    /// ```
    pub fn sequence_mut(&mut self) -> &mut Sequence {
        &mut self.sequence
    }

    /// Returns a mutable reference to the quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     self as sam,
    ///     record::{quality_scores, QualityScores},
    ///     AlignmentRecord,
    /// };
    ///
    /// let mut record = sam::Record::default();
    /// assert!(record.quality_scores().is_empty());
    ///
    /// let quality_scores: QualityScores = "NDLS".parse()?;
    /// *record.quality_scores_mut() = quality_scores.clone();
    ///
    /// assert_eq!(record.quality_scores(), &quality_scores);
    /// # Ok::<(), quality_scores::ParseError>(())
    /// ```
    pub fn quality_scores_mut(&mut self) -> &mut QualityScores {
        &mut self.quality_scores
    }

    /// Returns a mutable reference to the data fields for this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::data, AlignmentRecord};
    ///
    /// let mut record = sam::Record::default();
    /// assert!(record.data().is_empty());
    ///
    /// let field = data::Field::new(
    ///     data::field::Tag::AlignmentHitCount,
    ///     data::field::Value::Int32(1),
    /// );
    ///
    /// let data = record.data_mut();
    /// data.insert(field.clone());
    ///
    /// let data = record.data();
    /// assert_eq!(data.len(), 1);
    /// assert_eq!(data.get(field.tag()), Some(&field));
    /// ```
    pub fn data_mut(&mut self) -> &mut Data {
        &mut self.data
    }
}

impl AlignmentRecord for Record {
    fn read_name(&self) -> Option<&ReadName> {
        self.read_name.as_ref()
    }

    /// Returns the associated reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::ReferenceSequences, AlignmentRecord};
    /// let record = sam::Record::default();
    /// let reference_sequences = ReferenceSequences::default();
    /// assert!(record.reference_sequence(&reference_sequences).is_none());
    /// ```
    fn reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>> {
        get_reference_sequence(reference_sequences, self.reference_sequence_name())
    }

    fn flags(&self) -> Flags {
        self.flags
    }

    /// Returns the start position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, AlignmentRecord};
    /// let record = sam::Record::default();
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
    /// use noodles_sam::{self as sam, AlignmentRecord};
    /// let record = sam::Record::default();
    /// assert_eq!(record.alignment_span(), 0);
    /// ```
    fn alignment_span(&self) -> usize {
        self.cigar().reference_len()
    }

    fn mapping_quality(&self) -> Option<MappingQuality> {
        self.mapping_quality
    }

    fn cigar(&self) -> &Cigar {
        &self.cigar
    }

    /// Returns the associated reference sequence of the mate.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::ReferenceSequences, AlignmentRecord};
    /// let record = sam::Record::default();
    /// let reference_sequences = ReferenceSequences::default();
    /// assert!(record.mate_reference_sequence(&reference_sequences).is_none());
    /// ```
    fn mate_reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>> {
        get_reference_sequence(reference_sequences, self.mate_reference_sequence_name())
    }

    fn mate_alignment_start(&self) -> Option<Position> {
        self.mate_position()
    }

    fn template_length(&self) -> i32 {
        self.template_length
    }

    fn sequence(&self) -> Sequence {
        self.sequence.clone()
    }

    fn quality_scores(&self) -> &QualityScores {
        &self.quality_scores
    }

    fn data(&self) -> &Data {
        &self.data
    }
}

fn get_reference_sequence<'rs>(
    reference_sequences: &'rs ReferenceSequences,
    reference_sequence_name: Option<&ReferenceSequenceName>,
) -> Option<io::Result<&'rs ReferenceSequence>> {
    reference_sequence_name.map(|name| {
        reference_sequences.get(name.as_str()).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid reference sequence name",
            )
        })
    })
}

impl Default for Record {
    fn default() -> Self {
        Self {
            read_name: None,
            flags: Flags::UNMAPPED,
            reference_sequence_name: None,
            position: None,
            mapping_quality: None,
            cigar: Cigar::default(),
            mate_reference_sequence_name: None,
            mate_position: None,
            template_length: 0,
            sequence: Sequence::default(),
            quality_scores: QualityScores::default(),
            data: Data::default(),
        }
    }
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let qname = self
            .read_name()
            .map(|name| name.as_ref())
            .unwrap_or(NULL_FIELD);

        let rname = self
            .reference_sequence_name()
            .map(|name| name.as_str())
            .unwrap_or(NULL_FIELD);

        let pos = self.position().map(usize::from).unwrap_or_default();

        let mapq = self
            .mapping_quality()
            .map(u8::from)
            .unwrap_or(mapping_quality::MISSING);

        let rnext = self
            .mate_reference_sequence_name()
            .map(|mate_reference_sequence_name| {
                if let Some(reference_sequence_name) = self.reference_sequence_name() {
                    if mate_reference_sequence_name == reference_sequence_name {
                        return EQ_FIELD;
                    }
                }

                mate_reference_sequence_name.as_str()
            })
            .unwrap_or(NULL_FIELD);

        let pnext = self.mate_position().map(usize::from).unwrap_or_default();

        write!(
            f,
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}",
            qname = qname,
            flag = u16::from(self.flags()),
            rname = rname,
            pos = pos,
            mapq = mapq,
        )?;

        f.write_str("\t")?;
        if self.cigar().is_empty() {
            f.write_str(NULL_FIELD)?;
        } else {
            write!(f, "{}", self.cigar())?;
        }

        write!(
            f,
            "\t{rnext}\t{pnext}\t{tlen}",
            rnext = rnext,
            pnext = pnext,
            tlen = self.template_length(),
        )?;

        f.write_str("\t")?;
        if self.sequence().is_empty() {
            f.write_str(NULL_FIELD)?;
        } else {
            write!(f, "{}", self.sequence())?;
        }

        f.write_str("\t")?;
        if self.quality_scores().is_empty() {
            f.write_str(NULL_FIELD)?;
        } else {
            write!(f, "{}", self.quality_scores())?;
        }

        if !self.data().is_empty() {
            write!(f, "\t{}", self.data())?;
        }

        Ok(())
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parser::parse(s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let record = Record::default();

        assert!(record.read_name.is_none());
        assert_eq!(record.flags, Flags::UNMAPPED);
        assert!(record.reference_sequence_name.is_none());
        assert!(record.position.is_none());
        assert!(record.mapping_quality.is_none());
        assert!(record.cigar.is_empty());
        assert!(record.mate_reference_sequence_name.is_none());
        assert!(record.mate_position.is_none());
        assert_eq!(record.template_length, 0);
        assert!(record.sequence.is_empty());
        assert!(record.quality_scores.is_empty());
        assert!(record.data.is_empty());
    }

    #[test]
    fn test_fmt() {
        let record = Record::default();
        assert_eq!(record.to_string(), "*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*");
    }

    #[test]
    fn test_fmt_with_data() -> Result<(), Box<dyn std::error::Error>> {
        use super::data::{
            field::{Tag, Value},
            Field,
        };

        let data = Data::try_from(vec![Field::new(
            Tag::ReadGroup,
            Value::String(String::from("rg0")),
        )])?;

        let record = Record::builder().set_data(data).build();

        assert_eq!(
            record.to_string(),
            "*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\tRG:Z:rg0"
        );

        Ok(())
    }
}
