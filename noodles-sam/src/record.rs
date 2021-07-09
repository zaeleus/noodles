//! SAM record and fields.

pub mod builder;
pub mod cigar;
pub mod data;
mod field;
mod flags;
mod mapping_quality;
pub mod position;
pub mod quality_scores;
pub mod read_name;
pub mod reference_sequence_name;
pub mod sequence;

pub use self::{
    builder::Builder, cigar::Cigar, data::Data, field::Field, flags::Flags,
    mapping_quality::MappingQuality, position::Position, quality_scores::QualityScores,
    read_name::ReadName, reference_sequence_name::ReferenceSequenceName, sequence::Sequence,
};

use std::{error, fmt, num, str::FromStr};

pub(crate) const NULL_FIELD: &str = "*";
const ZERO_FIELD: &str = "0";
const EQ_FIELD: &str = "=";
const FIELD_DELIMITER: char = '\t';
const MAX_FIELDS: usize = 12;

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

    /// Returns the start position of this record.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
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

    /// Returns the start position of the mate of this record.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
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
    ///     .set_cigar("2M".parse()?)
    ///     .set_sequence("AT".parse()?)
    ///     .set_quality_scores("ND".parse()?)
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
    ///     .set_cigar("2M".parse()?)
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
    /// # use std::convert::TryFrom;
    /// use noodles_sam::{self as sam, record::{data, Data}};
    ///
    /// let record = sam::Record::default();
    /// assert!(record.data().is_empty());
    ///
    /// let data = Data::try_from(vec![data::Field::new(
    ///     data::field::Tag::AlignmentHitCount,
    ///     data::field::Value::Int32(1),
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
    ///     data::field::Value::Int32(1),
    /// );
    ///
    /// let data = record.data_mut();
    /// data.insert(field.tag().clone(), field.clone());
    ///
    /// let data = record.data();
    /// assert_eq!(data.len(), 1);
    /// assert_eq!(data.get(field.tag()), Some(&field));
    /// ```
    pub fn data_mut(&mut self) -> &mut Data {
        &mut self.data
    }
}

impl Default for Record {
    fn default() -> Self {
        // TODO
        Builder::new().build().unwrap()
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

/// An error returned when a raw SAM record fails to parse.
#[derive(Clone, Debug, PartialEq)]
pub enum ParseError {
    /// A required record field is missing.
    MissingField(Field),
    /// The record read name is invalid.
    InvalidReadName(read_name::ParseError),
    /// The record flags field is invalid.
    InvalidFlags(num::ParseIntError),
    /// The record reference sequence name is invalid.
    InvalidReferenceSequenceName(reference_sequence_name::ParseError),
    /// The record position is invalid.
    InvalidPosition(position::ParseError),
    /// The record mapping quality is invalid.
    InvalidMappingQuality(num::ParseIntError),
    /// The record CIGAR string is invalid.
    InvalidCigar(cigar::ParseError),
    /// The record mate reference sequence name is invalid.
    InvalidMateReferenceSequenceName(reference_sequence_name::ParseError),
    /// The record mate position is invalid.
    InvalidMatePosition(position::ParseError),
    /// The record template length is invalid.
    InvalidTemplateLength(num::ParseIntError),
    /// The record sequence is invalid.
    InvalidSequence(sequence::ParseError),
    /// The sequence length does not match the CIGAR string read length.
    SequenceLengthMismatch(u32, u32),
    /// The record quality score is invalid.
    InvalidQualityScores(quality_scores::ParseError),
    /// The quality scores length does not match the sequence length.
    QualityScoresLengthMismatch(u32, u32),
    /// The record data is invalid.
    InvalidData(data::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(field) => write!(f, "missing field: {}", field),
            Self::InvalidReadName(e) => write!(f, "invalid read name: {}", e),
            Self::InvalidFlags(e) => write!(f, "invalid flags: {}", e),
            Self::InvalidReferenceSequenceName(e) => {
                write!(f, "invalid reference sequence name: {}", e)
            }
            Self::InvalidPosition(e) => write!(f, "invalid position: {}", e),
            Self::InvalidMappingQuality(e) => write!(f, "invalid mapping quality: {}", e),
            Self::InvalidCigar(e) => write!(f, "invalid CIGAR: {}", e),
            Self::InvalidMateReferenceSequenceName(e) => {
                write!(f, "invalid mate reference sequence name: {}", e)
            }
            Self::InvalidMatePosition(e) => write!(f, "invalid mate position: {}", e),
            Self::InvalidTemplateLength(e) => write!(f, "invalid template length: {}", e),
            Self::InvalidSequence(e) => write!(f, "invalid sequence: {}", e),
            Self::SequenceLengthMismatch(sequence_len, cigar_read_len) => write!(
                f,
                "sequence length mismatch: expected {}, got {}",
                cigar_read_len, sequence_len
            ),
            Self::QualityScoresLengthMismatch(quality_scores_len, sequence_len) => write!(
                f,
                "quality scores length mismatch: expected {}, got {}",
                sequence_len, quality_scores_len
            ),
            Self::InvalidQualityScores(e) => write!(f, "invalid quality scores: {}", e),
            Self::InvalidData(e) => write!(f, "invalid data: {}", e),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use builder::BuildError;

        let mut fields = s.splitn(MAX_FIELDS, FIELD_DELIMITER);

        let mut builder = Self::builder();

        if let Some(qname) = parse_qname(&mut fields)? {
            builder = builder.set_read_name(qname);
        }

        let flag = parse_string(&mut fields, Field::Flags)
            .and_then(|s| s.parse::<u16>().map_err(ParseError::InvalidFlags))
            .map(Flags::from)?;

        builder = builder.set_flags(flag);

        let rname = parse_rname(&mut fields)?;

        if let Some(pos) = parse_pos(&mut fields)? {
            builder = builder.set_position(pos);
        }

        let mapq = parse_string(&mut fields, Field::MappingQuality)
            .and_then(|s| s.parse::<u8>().map_err(ParseError::InvalidMappingQuality))
            .map(MappingQuality::from)?;

        builder = builder.set_mapping_quality(mapq);

        let cigar: Cigar = parse_string(&mut fields, Field::Cigar)
            .and_then(|s| s.parse().map_err(ParseError::InvalidCigar))?;

        builder = builder.set_cigar(cigar);

        if let Some(rnext) = parse_rnext(&mut fields, rname.as_ref())? {
            builder = builder.set_mate_reference_sequence_name(rnext);
        }

        if let Some(reference_sequence_name) = rname {
            builder = builder.set_reference_sequence_name(reference_sequence_name);
        }

        if let Some(pnext) = parse_pnext(&mut fields)? {
            builder = builder.set_mate_position(pnext);
        }

        let tlen = parse_string(&mut fields, Field::TemplateLength)
            .and_then(|s| s.parse::<i32>().map_err(ParseError::InvalidTemplateLength))?;

        builder = builder.set_template_length(tlen);

        let seq = parse_string(&mut fields, Field::Sequence)
            .and_then(|s| s.parse().map_err(ParseError::InvalidSequence))?;

        builder = builder.set_sequence(seq);

        let qual = parse_string(&mut fields, Field::QualityScores)
            .and_then(|s| s.parse().map_err(ParseError::InvalidQualityScores))?;

        builder = builder.set_quality_scores(qual);

        if let Some(data) = parse_data(&mut fields)? {
            builder = builder.set_data(data);
        }

        match builder.build() {
            Ok(r) => Ok(r),
            Err(BuildError::SequenceLengthMismatch(sequence_len, cigar_read_len)) => Err(
                ParseError::SequenceLengthMismatch(sequence_len, cigar_read_len),
            ),
            Err(BuildError::QualityScoresLengthMismatch(quality_scores_len, sequence_len)) => Err(
                ParseError::QualityScoresLengthMismatch(quality_scores_len, sequence_len),
            ),
        }
    }
}

fn parse_string<'a, I>(fields: &mut I, field: Field) -> Result<&'a str, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields.next().ok_or(ParseError::MissingField(field))
}

fn parse_qname<'a, I>(fields: &mut I) -> Result<Option<ReadName>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    parse_string(fields, Field::Name).and_then(|s| {
        if s == NULL_FIELD {
            Ok(None)
        } else {
            s.parse().map(Some).map_err(ParseError::InvalidReadName)
        }
    })
}

fn parse_rname<'a, I>(fields: &mut I) -> Result<Option<ReferenceSequenceName>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    parse_string(fields, Field::ReferenceSequenceName).and_then(|s| {
        if s == NULL_FIELD {
            Ok(None)
        } else {
            s.parse()
                .map(Some)
                .map_err(ParseError::InvalidReferenceSequenceName)
        }
    })
}

fn parse_pos<'a, I>(fields: &mut I) -> Result<Option<Position>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    parse_string(fields, Field::Position).and_then(|s| match s {
        ZERO_FIELD => Ok(None),
        _ => s.parse().map(Some).map_err(ParseError::InvalidPosition),
    })
}

fn parse_rnext<'a, I>(
    fields: &mut I,
    rname: Option<&ReferenceSequenceName>,
) -> Result<Option<ReferenceSequenceName>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    parse_string(fields, Field::MateReferenceSequenceName).and_then(|s| match s {
        NULL_FIELD => Ok(None),
        EQ_FIELD => Ok(rname.cloned()),
        _ => s
            .parse()
            .map(Some)
            .map_err(ParseError::InvalidMateReferenceSequenceName),
    })
}

fn parse_pnext<'a, I>(fields: &mut I) -> Result<Option<Position>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    parse_string(fields, Field::MatePosition).and_then(|s| match s {
        ZERO_FIELD => Ok(None),
        _ => s.parse().map(Some).map_err(ParseError::InvalidMatePosition),
    })
}

fn parse_data<'a, I>(fields: &mut I) -> Result<Option<Data>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .map(|s| s.parse().map_err(ParseError::InvalidData))
        .transpose()
}

#[cfg(test)]
mod tests {
    use std::convert::TryFrom;

    use super::*;

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

    #[test]
    fn test_from_str_with_invalid_position() {
        let s = "*\t0\tsq0\t-1\t255\t4M\t*\t0\t0\tACGT\tNDLS";

        assert!(matches!(
            s.parse::<Record>(),
            Err(ParseError::InvalidPosition(_))
        ));

        let s = "*\t0\tsq0\tzero\t255\t4M\t*\t0\t0\tACGT\tNDLS";

        assert!(matches!(
            s.parse::<Record>(),
            Err(ParseError::InvalidPosition(_))
        ));
    }

    #[test]
    fn test_from_str_with_sequence_length_mismatch() {
        let s = "*\t0\tsq0\t1\t255\t2M\t*\t0\t0\tACGT\tNDLS";

        assert_eq!(
            s.parse::<Record>(),
            Err(ParseError::SequenceLengthMismatch(4, 2))
        );
    }

    #[test]
    fn test_from_str_with_quality_scores_length_mismatch() {
        let s = "*\t0\tsq0\t1\t255\t4M\t*\t0\t0\tACGT\tNDL";

        assert_eq!(
            s.parse::<Record>(),
            Err(ParseError::QualityScoresLengthMismatch(3, 4))
        );
    }
}
