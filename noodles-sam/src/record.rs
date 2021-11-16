//! SAM record and fields.

pub mod builder;
pub mod cigar;
pub mod data;
mod field;
mod flags;
mod mapping_quality;
mod parser;
pub mod position;
pub mod quality_scores;
pub mod read_name;
pub mod reference_sequence_name;
pub mod sequence;

pub use self::{
    builder::Builder, cigar::Cigar, data::Data, field::Field, flags::Flags,
    mapping_quality::MappingQuality, parser::ParseError, position::Position,
    quality_scores::QualityScores, read_name::ReadName,
    reference_sequence_name::ReferenceSequenceName, sequence::Sequence,
};

use std::{fmt, io, str::FromStr};

use super::{
    header::{ReferenceSequence, ReferenceSequences},
    RecordExt,
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
    mapping_quality: MappingQuality,
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
    /// use noodles_sam::{self as sam, record::Flags};
    ///
    /// let record = sam::Record::builder()
    ///     .set_read_name("r0".parse()?)
    ///     .set_flags(Flags::UNMAPPED)
    ///     .build()?;
    ///
    /// assert_eq!(record.read_name().map(|name| name.as_str()), Some("r0"));
    /// assert_eq!(record.flags(), Flags::UNMAPPED);
    /// assert!(record.reference_sequence_name().is_none());
    /// assert!(record.position().is_none());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn builder() -> Builder {
        Builder::new()
    }

    /// Returns the read name of this record.
    ///
    /// This is also called the query name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let record = sam::Record::default();
    /// assert!(record.read_name().is_none());
    ///
    /// let record = sam::Record::builder()
    ///     .set_read_name("r0".parse()?)
    ///     .build()?;
    /// assert_eq!(record.read_name().map(|name| name.as_str()), Some("r0"));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn read_name(&self) -> Option<&ReadName> {
        self.read_name.as_ref()
    }

    /// Returns a mutable reference to the read name.
    ///
    /// This is also called the query name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let mut record = sam::Record::default();
    /// record.read_name_mut().insert("r0".parse()?);
    /// assert_eq!(record.read_name().map(|name| name.as_str()), Some("r0"));
    ///
    /// *record.read_name_mut() = None;
    /// assert!(record.read_name().is_none());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn read_name_mut(&mut self) -> &mut Option<ReadName> {
        &mut self.read_name
    }

    /// Returns the SAM flags of this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Flags};
    ///
    /// let record = sam::Record::default();
    /// assert_eq!(record.flags(), Flags::UNMAPPED);
    /// assert_eq!(u16::from(record.flags()), 4);
    ///
    /// let record = sam::Record::builder()
    ///     .set_flags(Flags::PAIRED | Flags::READ_1)
    ///     .build()?;
    /// assert_eq!(record.flags(), Flags::PAIRED | Flags::READ_1);
    /// # Ok::<(), sam::record::builder::BuildError>(())
    /// ```
    pub fn flags(&self) -> Flags {
        self.flags
    }

    /// Returns a mutable reference to the SAM flags.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Flags};
    ///
    /// let mut record = sam::Record::builder()
    ///     .set_flags(Flags::PAIRED | Flags::READ_1)
    ///     .build()?;
    /// record.flags_mut().set(Flags::DUPLICATE, true);
    /// assert_eq!(record.flags(), Flags::PAIRED | Flags::READ_1 | Flags::DUPLICATE);
    ///
    /// record.flags_mut().set(Flags::PAIRED | Flags::QC_FAIL, false);
    /// assert_eq!(record.flags(), Flags::READ_1 | Flags::DUPLICATE);
    ///
    /// # Ok::<(), sam::record::builder::BuildError>(())
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
    /// assert_eq!(record.reference_sequence_name(), None);
    ///
    /// let record = sam::Record::builder()
    ///     .set_reference_sequence_name("sq0".parse()?)
    ///     .build()?;
    /// assert_eq!(record.reference_sequence_name().map(|name| name.as_str()), Some("sq0"));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
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
    /// use noodles_sam::{self as sam, record::Position};
    ///
    /// let record = sam::Record::default();
    /// assert!(record.position().is_none());
    ///
    /// let record = sam::Record::builder()
    ///     .set_position(Position::try_from(13)?)
    ///     .build()?;
    /// assert_eq!(record.position().map(i32::from), Some(13));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
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
    /// use noodles_sam::{self as sam, record::Position};
    ///
    /// let mut record = sam::Record::default();
    /// *record.position_mut() = Some(Position::try_from(13)?);
    /// assert_eq!(record.position().map(i32::from), Some(13));
    ///
    /// *record.position_mut() = None;
    /// assert!(record.position().is_none());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn position_mut(&mut self) -> &mut Option<Position> {
        &mut self.position
    }

    /// Returns the mapping quality of this record.
    ///
    /// Mapping quality ranges from 0 to 254, inclusive. A value of 255 means no mapping quality is
    /// set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::MappingQuality};
    ///
    /// let record = sam::Record::default();
    /// assert!(record.mapping_quality().is_none());
    /// assert_eq!(u8::from(record.mapping_quality()), 255);
    ///
    /// let record = sam::Record::builder().set_mapping_quality(MappingQuality::from(8)).build()?;
    /// assert_eq!(*record.mapping_quality(), Some(8));
    /// # Ok::<(), sam::record::builder::BuildError>(())
    /// ```
    pub fn mapping_quality(&self) -> MappingQuality {
        self.mapping_quality
    }

    /// Returns a mutable reference to the mapping quality.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::MappingQuality};
    ///
    /// let mut record = sam::Record::default();
    /// *record.mapping_quality_mut() = MappingQuality::from(8);
    /// assert_eq!(*record.mapping_quality(), Some(8));
    ///
    /// *record.mapping_quality_mut() = MappingQuality::from(255);
    /// assert!(record.mapping_quality().is_none());
    /// ```
    pub fn mapping_quality_mut(&mut self) -> &mut MappingQuality {
        &mut self.mapping_quality
    }

    /// Returns the CIGAR operations that describe how the read as mapped.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::cigar::{op, Op}};
    ///
    /// let record = sam::Record::default();
    /// assert!(record.cigar().is_empty());
    /// assert_eq!(record.cigar().to_string(), "*");
    ///
    /// let record = sam::Record::builder().set_cigar("34M2S".parse()?).build()?;
    /// assert_eq!(record.cigar().to_string(), "34M2S");
    ///
    /// let mut ops = record.cigar().iter();
    /// assert_eq!(ops.next(), Some(&Op::new(op::Kind::Match, 34)));
    /// assert_eq!(ops.next(), Some(&Op::new(op::Kind::SoftClip, 2)));
    /// assert_eq!(ops.next(), None);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn cigar(&self) -> &Cigar {
        &self.cigar
    }

    /// Returns a mutable reference to the CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::{cigar::{op, Op}, Cigar}};
    ///
    /// let mut record = sam::Record::default();
    /// assert!(record.cigar().is_empty());
    ///
    /// let cigar = Cigar::from(vec![
    ///     Op::new(op::Kind::Match, 36),
    ///     Op::new(op::Kind::SoftClip, 2),
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
    ///     .build()?;
    ///
    /// assert_eq!(
    ///     record.mate_reference_sequence_name().map(|name| name.as_str()),
    ///     Some("sq0")
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
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
    /// use noodles_sam::{self as sam, record::Position};
    ///
    /// let record = sam::Record::default();
    /// assert!(record.mate_position().is_none());
    ///
    /// let record = sam::Record::builder()
    ///     .set_mate_position(Position::try_from(21)?)
    ///     .build()?;
    /// assert_eq!(record.mate_position().map(i32::from), Some(21));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn mate_position(&self) -> Option<Position> {
        self.mate_position
    }

    /// Returns a mutable reference to the start position.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::Position};
    ///
    /// let mut record = sam::Record::default();
    /// *record.mate_position_mut() = Some(Position::try_from(13)?);
    /// assert_eq!(record.mate_position().map(i32::from), Some(13));
    ///
    /// *record.mate_position_mut() = None;
    /// assert!(record.mate_position().is_none());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn mate_position_mut(&mut self) -> &mut Option<Position> {
        &mut self.mate_position
    }

    /// Returns the template length of this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let record = sam::Record::default();
    /// assert_eq!(record.template_length(), 0);
    ///
    /// let record = sam::Record::builder().set_template_length(101).build()?;
    /// assert_eq!(record.template_length(), 101);
    /// # Ok::<(), sam::record::builder::BuildError>(())
    /// ```
    pub fn template_length(&self) -> i32 {
        self.template_length
    }

    /// Returns a mutable reference to the template length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::MappingQuality};
    ///
    /// let mut record = sam::Record::default();
    /// *record.template_length_mut() = 101;
    /// assert_eq!(record.template_length(), 101);
    /// ```
    pub fn template_length_mut(&mut self) -> &mut i32 {
        &mut self.template_length
    }

    /// Returns the bases in the sequence of this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::sequence::Base};
    ///
    /// let record = sam::Record::default();
    /// assert!(record.sequence().is_empty());
    /// assert_eq!(record.sequence().to_string(), "*");
    ///
    /// let record = sam::Record::builder()
    ///     .set_sequence("AT".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(record.sequence().to_string(), "AT");
    ///
    /// let mut bases = record.sequence().iter();
    /// assert_eq!(bases.next(), Some(&Base::A));
    /// assert_eq!(bases.next(), Some(&Base::T));
    /// assert_eq!(bases.next(), None);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn sequence(&self) -> &Sequence {
        &self.sequence
    }

    /// Returns a mutable reference to the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::{sequence, Sequence}};
    ///
    /// let mut record = sam::Record::default();
    /// assert!(record.sequence().is_empty());
    ///
    /// let sequence: Sequence = "ACGT".parse()?;
    /// *record.sequence_mut() = sequence.clone();
    ///
    /// assert_eq!(record.sequence(), &sequence);
    /// # Ok::<(), sequence::ParseError>(())
    /// ```
    pub fn sequence_mut(&mut self) -> &mut Sequence {
        &mut self.sequence
    }

    /// Returns the quality score for each base in the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::quality_scores::Score};
    ///
    /// let record = sam::Record::default();
    /// assert!(record.quality_scores().is_empty());
    /// assert_eq!(record.quality_scores().to_string(), "*");
    ///
    /// let record = sam::Record::builder()
    ///     .set_sequence("AC".parse()?)
    ///     .set_quality_scores("ND".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(record.quality_scores().to_string(), "ND");
    ///
    /// let mut scores = record.quality_scores().iter().copied().map(u8::from);
    /// assert_eq!(scores.next(), Some(45));
    /// assert_eq!(scores.next(), Some(35));
    /// assert_eq!(scores.next(), None);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn quality_scores(&self) -> &QualityScores {
        &self.quality_scores
    }

    /// Returns a mutable reference to the quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::{quality_scores, QualityScores}};
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

    /// Returns the optional data fields for this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::{data, Data}};
    ///
    /// let record = sam::Record::default();
    /// assert!(record.data().is_empty());
    ///
    /// let data = Data::try_from(vec![data::Field::new(
    ///     data::field::Tag::AlignmentHitCount,
    ///     data::field::Value::Int(1),
    /// )])?;
    /// let record = sam::Record::builder().set_data(data).build()?;
    /// assert_eq!(record.data().to_string(), "NH:i:1");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn data(&self) -> &Data {
        &self.data
    }

    /// Returns a mutable reference to the data fields for this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, record::data};
    ///
    /// let mut record = sam::Record::default();
    /// assert!(record.data().is_empty());
    ///
    /// let field = data::Field::new(
    ///     data::field::Tag::AlignmentHitCount,
    ///     data::field::Value::Int(1),
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

impl RecordExt for Record {
    /// Returns the reference sequence name from the referece sequence dictionary.
    fn reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>> {
        self.reference_sequence_name()
            .map(|reference_sequence_name| {
                reference_sequences
                    .get(reference_sequence_name.as_str())
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "invalid reference sequence name",
                        )
                    })
            })
    }
}

impl Default for Record {
    fn default() -> Self {
        Self {
            read_name: Default::default(),
            flags: Flags::UNMAPPED,
            reference_sequence_name: Default::default(),
            position: Default::default(),
            mapping_quality: MappingQuality::default(),
            cigar: Cigar::default(),
            mate_reference_sequence_name: Default::default(),
            mate_position: Default::default(),
            template_length: Default::default(),
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
            .map(|name| name.as_str())
            .unwrap_or(NULL_FIELD);

        let rname = self
            .reference_sequence_name()
            .map(|name| name.as_str())
            .unwrap_or(NULL_FIELD);

        let pos = self.position().map(i32::from).unwrap_or(position::UNMAPPED);

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

        let pnext = self
            .mate_position()
            .map(i32::from)
            .unwrap_or(position::UNMAPPED);

        write!(
            f,
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}",
            qname = qname,
            flag = u16::from(self.flags()),
            rname = rname,
            pos = pos,
            mapq = u8::from(self.mapping_quality()),
            cigar = self.cigar(),
            rnext = rnext,
            pnext = pnext,
            tlen = self.template_length(),
            seq = self.sequence(),
            qual = self.quality_scores(),
        )?;

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
        let data = Data::try_from(vec![data::Field::new(
            data::field::Tag::ReadGroup,
            data::field::Value::String(String::from("rg0")),
        )])?;

        let record = Record::builder().set_data(data).build()?;

        assert_eq!(
            record.to_string(),
            "*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\tRG:Z:rg0"
        );

        Ok(())
    }
}
