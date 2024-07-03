//! BED record and fields.

pub mod builder;
pub mod color;
pub mod name;
pub mod score;
pub mod strand;

pub use self::{builder::Builder, color::Color, name::Name, score::Score, strand::Strand};

use std::{
    error,
    fmt::{self, Write},
    num::{self, NonZeroUsize},
    ops::Deref,
    str::FromStr,
};

use noodles_core::Position;

const DELIMITER: char = '\t';
const MISSING_STRING: &str = ".";
const MISSING_NUMBER: &str = "0";

type Block = (usize, usize);

#[derive(Clone, Debug, Eq, PartialEq)]
struct StandardFields {
    reference_sequence_name: String,
    start_position: Position,
    end_position: Position,
    name: Option<Name>,
    score: Option<Score>,
    strand: Option<Strand>,
    thick_start: Position,
    thick_end: Position,
    color: Option<Color>,
    blocks: Vec<Block>,
}

impl StandardFields {
    fn new<N>(reference_sequence_name: N, start_position: Position, end_position: Position) -> Self
    where
        N: Into<String>,
    {
        Self {
            reference_sequence_name: reference_sequence_name.into(),
            start_position,
            end_position,
            name: None,
            score: None,
            strand: None,
            thick_start: start_position,
            thick_end: end_position,
            color: None,
            blocks: Vec::new(),
        }
    }
}

/// Raw BED record optional fields.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct OptionalFields(Vec<String>);

impl Deref for OptionalFields {
    type Target = [String];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for OptionalFields {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, field) in self.0.iter().enumerate() {
            if i > 0 {
                f.write_char(DELIMITER)?;
            }

            f.write_str(field)?;
        }

        Ok(())
    }
}

impl From<Vec<String>> for OptionalFields {
    fn from(fields: Vec<String>) -> Self {
        Self(fields)
    }
}

#[cfg(test)]
mod optional_fields_tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let fields = OptionalFields::default();
        assert_eq!(fields.to_string(), "");

        let fields = OptionalFields::from(vec![String::from("n")]);
        assert_eq!(fields.to_string(), "n");

        let fields = OptionalFields::from(vec![String::from("n"), String::from("d")]);
        assert_eq!(fields.to_string(), "n\td");
    }
}

/// A BED record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RecordBuf<const N: u8> {
    standard_fields: StandardFields,
    optional_fields: OptionalFields,
}

/// A trait describing the number of standard fields in a BED record.
pub trait BedN<const N: u8> {}

impl BedN<3> for RecordBuf<3> {}
impl BedN<3> for RecordBuf<4> {}
impl BedN<3> for RecordBuf<5> {}
impl BedN<3> for RecordBuf<6> {}
impl BedN<3> for RecordBuf<7> {}
impl BedN<3> for RecordBuf<8> {}
impl BedN<3> for RecordBuf<9> {}
impl BedN<3> for RecordBuf<12> {}

impl BedN<4> for RecordBuf<4> {}
impl BedN<4> for RecordBuf<5> {}
impl BedN<4> for RecordBuf<6> {}
impl BedN<4> for RecordBuf<7> {}
impl BedN<4> for RecordBuf<8> {}
impl BedN<4> for RecordBuf<9> {}
impl BedN<4> for RecordBuf<12> {}

impl BedN<5> for RecordBuf<5> {}
impl BedN<5> for RecordBuf<6> {}
impl BedN<5> for RecordBuf<7> {}
impl BedN<5> for RecordBuf<8> {}
impl BedN<5> for RecordBuf<9> {}
impl BedN<5> for RecordBuf<12> {}

impl BedN<6> for RecordBuf<6> {}
impl BedN<6> for RecordBuf<7> {}
impl BedN<6> for RecordBuf<8> {}
impl BedN<6> for RecordBuf<9> {}
impl BedN<6> for RecordBuf<12> {}

impl BedN<7> for RecordBuf<7> {}
impl BedN<7> for RecordBuf<8> {}
impl BedN<7> for RecordBuf<9> {}
impl BedN<7> for RecordBuf<12> {}

impl BedN<8> for RecordBuf<8> {}
impl BedN<8> for RecordBuf<9> {}
impl BedN<8> for RecordBuf<12> {}

impl BedN<9> for RecordBuf<9> {}
impl BedN<9> for RecordBuf<12> {}

impl BedN<12> for RecordBuf<12> {}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<3>,
{
    /// Creates a BED record builder.
    pub fn builder() -> Builder<N> {
        Builder::default()
    }

    fn new(standard_fields: StandardFields, optional_fields: OptionalFields) -> Self {
        Self {
            standard_fields,
            optional_fields,
        }
    }

    /// Returns the reference sequence name (`chrom`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<3>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.reference_sequence_name(), "sq0");
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_sequence_name(&self) -> &str {
        &self.standard_fields.reference_sequence_name
    }

    /// Returns the feature start position (`chromStart`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let start_position = Position::try_from(8)?;
    ///
    /// let record = bed::feature::RecordBuf::<3>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(start_position)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.start_position(), start_position);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn start_position(&self) -> Position {
        self.standard_fields.start_position
    }

    /// Returns the feature end position (`chromEnd`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let end_position = Position::try_from(13)?;
    ///
    /// let record = bed::feature::RecordBuf::<3>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(end_position)
    ///     .build()?;
    ///
    /// assert_eq!(record.end_position(), end_position);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn end_position(&self) -> Position {
        self.standard_fields.end_position
    }

    /// Returns the list of raw optional fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<3>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .build()?;
    ///
    /// assert!(record.optional_fields().is_empty());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn optional_fields(&self) -> &OptionalFields {
        &self.optional_fields
    }
}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<4>,
{
    /// Returns the feature name (`name`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::{self as bed, feature::record_buf::Name};
    /// use noodles_core::Position;
    ///
    /// let name = Name::from("ndls1");
    ///
    /// let record = bed::feature::RecordBuf::<4>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_name(name.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.name(), Some(&name));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn name(&self) -> Option<&Name> {
        self.standard_fields.name.as_ref()
    }
}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<5>,
{
    /// Returns the score (`score`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::{self as bed, feature::record_buf::Score};
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<5>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_score(Score::try_from(21)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.score().map(u16::from), Some(21));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn score(&self) -> Option<Score> {
        self.standard_fields.score
    }
}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<6>,
{
    /// Returns the feature strand (`strand`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::{self as bed, feature::record_buf::Strand};
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<6>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_strand(Strand::Forward)
    ///     .build()?;
    ///
    /// assert_eq!(record.strand(), Some(Strand::Forward));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn strand(&self) -> Option<Strand> {
        self.standard_fields.strand
    }
}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<7>,
{
    /// Returns the thick start position (`thickStart`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let thick_start = Position::try_from(8)?;
    ///
    /// let record = bed::feature::RecordBuf::<7>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_thick_start(thick_start)
    ///     .build()?;
    ///
    /// assert_eq!(record.thick_start(), thick_start);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn thick_start(&self) -> Position {
        self.standard_fields.thick_start
    }
}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<8>,
{
    /// Returns the thick end position (`thickEnd`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let thick_end = Position::try_from(13)?;
    ///
    /// let record = bed::feature::RecordBuf::<8>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_thick_end(thick_end)
    ///     .build()?;
    ///
    /// assert_eq!(record.thick_end(), thick_end);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn thick_end(&self) -> Position {
        self.standard_fields.thick_end
    }
}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<9>,
{
    /// Returns the color (`itemRgb`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::{self as bed, feature::record_buf::Color};
    /// use noodles_core::Position;
    ///
    /// let thick_end = Position::try_from(13)?;
    ///
    /// let record = bed::feature::RecordBuf::<9>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_color(Color::RED)
    ///     .build()?;
    ///
    /// assert_eq!(record.color(), Some(Color::RED));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn color(&self) -> Option<Color> {
        self.standard_fields.color
    }
}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<12>,
{
    /// Returns the blocks (`[(blockStarts, blockSizes)]`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let blocks = vec![(0, 2)];
    ///
    /// let record = bed::feature::RecordBuf::<12>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_blocks(blocks.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.blocks(), &blocks);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn blocks(&self) -> &[Block] {
        &self.standard_fields.blocks
    }
}

/// An error returned when a raw BED record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The reference sequence name is missing.
    MissingReferenceSequenceName,
    /// The start position is missing.
    MissingStartPosition,
    /// The start position is invalid.
    InvalidStartPosition,
    /// The end position is missing.
    MissingEndPosition,
    /// The end position is invalid.
    InvalidEndPosition(num::ParseIntError),
    /// The name is missing.
    MissingName,
    /// The score is missing.
    MissingScore,
    /// The score is invalid.
    InvalidScore(score::ParseError),
    /// The strand is missing.
    MissingStrand,
    /// The strand is invalid.
    InvalidStrand(strand::ParseError),
    /// The thick start position is missing.
    MissingThickStart,
    /// The thick start position is invalid.
    InvalidThickStart,
    /// The thick end position is missing.
    MissingThickEnd,
    /// The thick end position is invalid.
    InvalidThickEnd(num::ParseIntError),
    /// The color is missing.
    MissingColor,
    /// The color is invalid.
    InvalidColor(color::ParseError),
    /// The block count is missing.
    MissingBlockCount,
    /// The block count is invalid.
    InvalidBlockCount(num::ParseIntError),
    /// The block sizes are missing.
    MissingBlockSizes,
    /// A block size is invalid.
    InvalidBlockSize(num::ParseIntError),
    /// The block starts are missing.
    MissingBlockStarts,
    /// A block start is invalid.
    InvalidBlockStart(num::ParseIntError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidEndPosition(e) | Self::InvalidThickEnd(e) | Self::InvalidBlockStart(e) => {
                Some(e)
            }
            Self::InvalidScore(e) => Some(e),
            Self::InvalidStrand(e) => Some(e),
            Self::InvalidColor(e) => Some(e),
            Self::InvalidBlockCount(e) | Self::InvalidBlockSize(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingReferenceSequenceName => f.write_str("missing reference sequence name"),
            Self::MissingStartPosition => f.write_str("missing start position"),
            Self::InvalidStartPosition => f.write_str("invalid start position"),
            Self::MissingEndPosition => f.write_str("missing end position"),
            Self::InvalidEndPosition(_) => f.write_str("invalid end position"),
            Self::MissingName => f.write_str("missing name"),
            Self::MissingScore => f.write_str("missing score"),
            Self::InvalidScore(_) => f.write_str("invalid score"),
            Self::MissingStrand => f.write_str("missing strand"),
            Self::InvalidStrand(_) => f.write_str("invalid strand"),
            Self::MissingThickStart => f.write_str("missing thick start"),
            Self::InvalidThickStart => f.write_str("invalid thick start"),
            Self::MissingThickEnd => f.write_str("missing thick end"),
            Self::InvalidThickEnd(_) => f.write_str("invalid thick end"),
            Self::MissingColor => f.write_str("missing color"),
            Self::InvalidColor(_) => f.write_str("invalid color"),
            Self::MissingBlockCount => f.write_str("missing block count"),
            Self::InvalidBlockCount(_) => f.write_str("invalid block count"),
            Self::MissingBlockSizes => f.write_str("missing block sizes"),
            Self::InvalidBlockSize(_) => f.write_str("invalid block size"),
            Self::MissingBlockStarts => f.write_str("missing block starts"),
            Self::InvalidBlockStart(_) => f.write_str("invalid block start"),
        }
    }
}

impl FromStr for RecordBuf<3> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        let standard_fields = parse_bed_3_fields(&mut fields)?;
        let optional_fields = parse_optional_fields(&mut fields);
        Ok(Self::new(standard_fields, optional_fields))
    }
}

impl FromStr for RecordBuf<4> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        let standard_fields = parse_bed_4_fields(&mut fields)?;
        let optional_fields = parse_optional_fields(&mut fields);
        Ok(Self::new(standard_fields, optional_fields))
    }
}

impl FromStr for RecordBuf<5> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        let standard_fields = parse_bed_5_fields(&mut fields)?;
        let optional_fields = parse_optional_fields(&mut fields);
        Ok(Self::new(standard_fields, optional_fields))
    }
}

impl FromStr for RecordBuf<6> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        let standard_fields = parse_bed_6_fields(&mut fields)?;
        let optional_fields = parse_optional_fields(&mut fields);
        Ok(Self::new(standard_fields, optional_fields))
    }
}

impl FromStr for RecordBuf<7> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        let standard_fields = parse_bed_7_fields(&mut fields)?;
        let optional_fields = parse_optional_fields(&mut fields);
        Ok(Self::new(standard_fields, optional_fields))
    }
}

impl FromStr for RecordBuf<8> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        let standard_fields = parse_bed_8_fields(&mut fields)?;
        let optional_fields = parse_optional_fields(&mut fields);
        Ok(Self::new(standard_fields, optional_fields))
    }
}

impl FromStr for RecordBuf<9> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        let standard_fields = parse_bed_9_fields(&mut fields)?;
        let optional_fields = parse_optional_fields(&mut fields);
        Ok(Self::new(standard_fields, optional_fields))
    }
}

impl FromStr for RecordBuf<12> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        let standard_fields = parse_bed_12_fields(&mut fields)?;
        let optional_fields = parse_optional_fields(&mut fields);
        Ok(Self::new(standard_fields, optional_fields))
    }
}

fn parse_bed_3_fields<'a, I>(fields: &mut I) -> Result<StandardFields, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    parse_mandatory_fields(fields)
}

fn parse_bed_4_fields<'a, I>(fields: &mut I) -> Result<StandardFields, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    let mut standard_fields = parse_bed_3_fields(fields)?;
    standard_fields.name = parse_name(fields)?;
    Ok(standard_fields)
}

fn parse_bed_5_fields<'a, I>(fields: &mut I) -> Result<StandardFields, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    let mut standard_fields = parse_bed_4_fields(fields)?;
    standard_fields.score = parse_score(fields)?;
    Ok(standard_fields)
}

fn parse_bed_6_fields<'a, I>(fields: &mut I) -> Result<StandardFields, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    let mut standard_fields = parse_bed_5_fields(fields)?;
    standard_fields.strand = parse_strand(fields)?;
    Ok(standard_fields)
}

fn parse_bed_7_fields<'a, I>(fields: &mut I) -> Result<StandardFields, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    let mut standard_fields = parse_bed_6_fields(fields)?;
    standard_fields.thick_start = parse_thick_start(fields)?;
    Ok(standard_fields)
}

fn parse_bed_8_fields<'a, I>(fields: &mut I) -> Result<StandardFields, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    let mut standard_fields = parse_bed_7_fields(fields)?;
    standard_fields.thick_end = parse_thick_end(fields)?;
    Ok(standard_fields)
}

fn parse_bed_9_fields<'a, I>(fields: &mut I) -> Result<StandardFields, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    let mut standard_fields = parse_bed_8_fields(fields)?;
    standard_fields.color = parse_color(fields)?;
    Ok(standard_fields)
}

fn parse_bed_12_fields<'a, I>(fields: &mut I) -> Result<StandardFields, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    let mut standard_fields = parse_bed_9_fields(fields)?;
    standard_fields.blocks = parse_blocks(fields)?;
    Ok(standard_fields)
}

fn parse_mandatory_fields<'a, I>(fields: &mut I) -> Result<StandardFields, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    let reference_sequence_name = fields
        .next()
        .ok_or(ParseError::MissingReferenceSequenceName)?;

    let start_position = fields
        .next()
        .ok_or(ParseError::MissingStartPosition)
        .and_then(|s| {
            s.parse()
                .map_err(|_| ParseError::InvalidStartPosition)
                .and_then(|n: usize| {
                    n.checked_add(1)
                        .ok_or(ParseError::InvalidStartPosition)
                        .and_then(|m| {
                            Position::try_from(m).map_err(|_| ParseError::InvalidStartPosition)
                        })
                })
        })?;

    let end_position = fields
        .next()
        .ok_or(ParseError::MissingEndPosition)
        .and_then(|s| s.parse().map_err(ParseError::InvalidEndPosition))?;

    Ok(StandardFields::new(
        reference_sequence_name,
        start_position,
        end_position,
    ))
}

fn parse_name<'a, I>(fields: &mut I) -> Result<Option<Name>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .map(|s| match s {
            MISSING_STRING => None,
            _ => Some(Name::from(s)),
        })
        .ok_or(ParseError::MissingName)
}

fn parse_score<'a, I>(fields: &mut I) -> Result<Option<Score>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingScore)
        .and_then(|s| match s {
            MISSING_NUMBER => Ok(None),
            _ => s.parse().map(Some).map_err(ParseError::InvalidScore),
        })
}

fn parse_strand<'a, I>(fields: &mut I) -> Result<Option<Strand>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingStrand)
        .and_then(|s| match s {
            MISSING_STRING => Ok(None),
            _ => s.parse().map(Some).map_err(ParseError::InvalidStrand),
        })
}

fn parse_thick_start<'a, I>(fields: &mut I) -> Result<Position, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingThickStart)
        .and_then(|s| {
            s.parse()
                .map_err(|_| ParseError::InvalidThickStart)
                .and_then(|n: usize| {
                    n.checked_add(1)
                        .ok_or(ParseError::InvalidThickStart)
                        .and_then(|m| {
                            Position::try_from(m).map_err(|_| ParseError::InvalidThickStart)
                        })
                })
        })
}

fn parse_thick_end<'a, I>(fields: &mut I) -> Result<Position, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingThickEnd)
        .and_then(|s| s.parse().map_err(ParseError::InvalidThickEnd))
}

fn parse_color<'a, I>(fields: &mut I) -> Result<Option<Color>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingColor)
        .and_then(|s| match s {
            MISSING_NUMBER => Ok(None),
            _ => s.parse().map(Some).map_err(ParseError::InvalidColor),
        })
}

fn parse_blocks<'a, I>(fields: &mut I) -> Result<Vec<Block>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    const LIST_DELIMITER: char = ',';

    let len = fields
        .next()
        .ok_or(ParseError::MissingBlockCount)
        .and_then(|s| {
            s.parse::<NonZeroUsize>()
                .map_err(ParseError::InvalidBlockCount)
        })
        .map(usize::from)?;

    let raw_sizes = fields
        .next()
        .ok_or(ParseError::MissingBlockSizes)
        .map(|s| s.splitn(len, LIST_DELIMITER))?;

    let raw_starts = fields
        .next()
        .ok_or(ParseError::MissingBlockStarts)
        .map(|s| s.splitn(len, LIST_DELIMITER))?;

    let mut blocks = Vec::with_capacity(len);

    for (raw_start, raw_size) in raw_starts.zip(raw_sizes) {
        let start = raw_start.parse().map_err(ParseError::InvalidBlockStart)?;
        let size = raw_size.parse().map_err(ParseError::InvalidBlockSize)?;
        blocks.push((start, size));
    }

    Ok(blocks)
}

fn parse_optional_fields<'a, I>(fields: &mut I) -> OptionalFields
where
    I: Iterator<Item = &'a str>,
{
    let raw_fields: Vec<_> = fields.map(|s| s.into()).collect();
    OptionalFields::from(raw_fields)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str_for_record_3() -> Result<(), noodles_core::position::TryFromIntError> {
        let actual = "sq0\t7\t13".parse::<RecordBuf<3>>();

        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        let standard_fields = StandardFields::new("sq0", start, end);
        let expected = Ok(RecordBuf::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_from_str_for_record_4() -> Result<(), Box<dyn std::error::Error>> {
        let actual = "sq0\t7\t13\tndls1".parse::<RecordBuf<4>>();

        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        let mut standard_fields = StandardFields::new("sq0", start, end);
        standard_fields.name = Some(Name::from("ndls1"));

        let expected = Ok(RecordBuf::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_from_str_for_record_5() -> Result<(), Box<dyn std::error::Error>> {
        let actual = "sq0\t7\t13\t.\t21".parse::<RecordBuf<5>>();

        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        let mut standard_fields = StandardFields::new("sq0", start, end);
        standard_fields.score = Score::try_from(21).map(Some)?;

        let expected = Ok(RecordBuf::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_from_str_for_record_6() -> Result<(), noodles_core::position::TryFromIntError> {
        let actual = "sq0\t7\t13\t.\t0\t+".parse::<RecordBuf<6>>();

        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        let mut standard_fields = StandardFields::new("sq0", start, end);
        standard_fields.strand = Some(Strand::Forward);

        let expected = Ok(RecordBuf::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_from_str_for_record_7() -> Result<(), noodles_core::position::TryFromIntError> {
        let actual = "sq0\t7\t13\t.\t0\t.\t7".parse::<RecordBuf<7>>();

        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        let mut standard_fields = StandardFields::new("sq0", start, end);
        standard_fields.thick_start = start;

        let expected = Ok(RecordBuf::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_from_str_for_record_8() -> Result<(), noodles_core::position::TryFromIntError> {
        let actual = "sq0\t7\t13\t.\t0\t.\t7\t13".parse::<RecordBuf<8>>();

        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        let mut standard_fields = StandardFields::new("sq0", start, end);
        standard_fields.thick_start = start;
        standard_fields.thick_end = end;

        let expected = Ok(RecordBuf::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_from_str_for_record_9() -> Result<(), noodles_core::position::TryFromIntError> {
        let actual = "sq0\t7\t13\t.\t0\t.\t7\t13\t255,0,0".parse::<RecordBuf<9>>();

        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        let mut standard_fields = StandardFields::new("sq0", start, end);
        standard_fields.thick_start = start;
        standard_fields.thick_end = end;
        standard_fields.color = Some(Color::RED);

        let expected = Ok(RecordBuf::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_from_str_for_record_12() -> Result<(), noodles_core::position::TryFromIntError> {
        let actual = "sq0\t7\t13\t.\t0\t.\t7\t13\t0\t1\t2\t0".parse::<RecordBuf<12>>();

        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        let mut standard_fields = StandardFields::new("sq0", start, end);
        standard_fields.thick_start = start;
        standard_fields.thick_end = end;
        standard_fields.blocks = vec![(0, 2)];

        let expected = Ok(RecordBuf::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);

        Ok(())
    }
}
