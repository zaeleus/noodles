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
    num,
    ops::Deref,
    str::FromStr,
};
use std::cmp::min;

use noodles_core::Position;

const DELIMITER: char = '\t';
const MISSING_STRING: &str = ".";
const MISSING_NUMBER: &str = "0";
const COMMA_SEPARATOR: char = ',';

#[derive(Clone, Debug, Eq, PartialEq)]
struct StandardFields {
    reference_sequence_name: String,
    start_position: Position,
    end_position: Position,
    name: Option<Name>,
    score: Option<Score>,
    strand: Option<Strand>,
    thick_start: Option<Position>,
    thick_end: Option<Position>,
    item_rgb: Option<Color>,
    block_sizes: Option<Vec<usize>>,
    block_starts: Option<Vec<Position>>,
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
            thick_start: None,
            thick_end: None,
            item_rgb: None,
            block_sizes: None,
            block_starts: None,
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
pub struct Record<const N: u8> {
    standard_fields: StandardFields,
    optional_fields: OptionalFields,
}

/// A trait describing the number of standard fields in a BED record.
pub trait BedN<const N: u8> {}

impl BedN<3> for Record<3> {}
impl BedN<3> for Record<4> {}
impl BedN<3> for Record<5> {}
impl BedN<3> for Record<6> {}
impl BedN<3> for Record<12> {}

impl BedN<4> for Record<4> {}
impl BedN<4> for Record<5> {}
impl BedN<4> for Record<6> {}
impl BedN<4> for Record<12> {}

impl BedN<5> for Record<5> {}
impl BedN<5> for Record<6> {}
impl BedN<5> for Record<12> {}

impl BedN<6> for Record<6> {}
impl BedN<6> for Record<12> {}

impl BedN<12> for Record<12> {}

impl<const N: u8> Record<N>
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
    /// let record = bed::Record::<3>::builder()
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
    /// let record = bed::Record::<3>::builder()
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
    /// let record = bed::Record::<3>::builder()
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
    /// let record = bed::Record::<3>::builder()
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

impl<const N: u8> Record<N>
where
    Self: BedN<4>,
{
    /// Returns the feature name (`name`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::{self as bed, record::Name};
    /// use noodles_core::Position;
    ///
    /// let name: Name = "ndls1".parse()?;
    ///
    /// let record = bed::Record::<4>::builder()
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

impl<const N: u8> Record<N>
where
    Self: BedN<5>,
{
    /// Returns the score (`score`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::{self as bed, record::Score};
    /// use noodles_core::Position;
    ///
    /// let record = bed::Record::<5>::builder()
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

impl<const N: u8> Record<N>
where
    Self: BedN<6>,
{
    /// Returns the feature strand (`strand`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::{self as bed, record::Strand};
    /// use noodles_core::Position;
    ///
    /// let record = bed::Record::<6>::builder()
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

impl<const N: u8> Record<N>
where
    Self: BedN<12>,
{
    pub fn thick_start(&self) -> Option<Position> {
        self.standard_fields.thick_start
    }

    pub fn think_end(&self) -> Option<Position> {
        self.standard_fields.thick_end
    }

    pub fn item_rgb(&self) -> Option<&Color> {
        self.standard_fields.item_rgb.as_ref()
    }

    pub fn block_sizes(&self) -> Option<&[usize]> {
        self.standard_fields
            .block_sizes
            .as_ref()
            .map(|i| i.as_slice())
    }

    pub fn block_starts(&self) -> Option<&[Position]> {
        self.standard_fields
            .block_starts
            .as_ref()
            .map(|i| i.as_slice())
    }
}

impl fmt::Display for Record<3> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        format_bed_3_fields(f, self)?;
        format_optional_fields(f, self.optional_fields())?;
        Ok(())
    }
}

impl fmt::Display for Record<4> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        format_bed_4_fields(f, self)?;
        format_optional_fields(f, self.optional_fields())?;
        Ok(())
    }
}

impl fmt::Display for Record<5> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        format_bed_5_fields(f, self)?;
        format_optional_fields(f, self.optional_fields())?;
        Ok(())
    }
}

impl fmt::Display for Record<6> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        format_bed_6_fields(f, self)?;
        format_optional_fields(f, self.optional_fields())?;
        Ok(())
    }
}

fn format_bed_3_fields<const N: u8>(f: &mut fmt::Formatter<'_>, record: &Record<N>) -> fmt::Result
where
    Record<N>: BedN<3>,
{
    write!(
        f,
        "{}{}{}{}{}",
        record.reference_sequence_name(),
        DELIMITER,
        usize::from(record.start_position()) - 1,
        DELIMITER,
        record.end_position()
    )
}

fn format_bed_4_fields<const N: u8>(f: &mut fmt::Formatter<'_>, record: &Record<N>) -> fmt::Result
where
    Record<N>: BedN<3> + BedN<4>,
{
    format_bed_3_fields(f, record)?;

    f.write_char(DELIMITER)?;

    if let Some(name) = record.name() {
        write!(f, "{}", name)
    } else {
        f.write_str(MISSING_STRING)
    }
}

fn format_bed_5_fields<const N: u8>(f: &mut fmt::Formatter<'_>, record: &Record<N>) -> fmt::Result
where
    Record<N>: BedN<3> + BedN<4> + BedN<5>,
{
    format_bed_4_fields(f, record)?;

    f.write_char(DELIMITER)?;

    if let Some(score) = record.score() {
        write!(f, "{}", score)
    } else {
        f.write_str(MISSING_NUMBER)
    }
}

fn format_bed_6_fields<const N: u8>(f: &mut fmt::Formatter<'_>, record: &Record<N>) -> fmt::Result
where
    Record<N>: BedN<3> + BedN<4> + BedN<5> + BedN<6>,
{
    format_bed_5_fields(f, record)?;

    f.write_char(DELIMITER)?;

    if let Some(strand) = record.strand() {
        write!(f, "{}", strand)
    } else {
        f.write_str(MISSING_STRING)
    }
}

fn format_bed_12_fields(f: &mut fmt::Formatter<'_>, record: &Record<12>) -> fmt::Result {
    format_bed_6_fields(f, record)?;

    f.write_char(DELIMITER)?;
    if let Some(think_start) = record.thick_start() {
        write!(f, "{}", think_start)?;
    } else {
        f.write_str(MISSING_STRING)?;
    }

    f.write_char(DELIMITER)?;
    if let Some(thick_end) = record.think_end() {
        write!(f, "{}", thick_end)?;
    } else {
        f.write_str(MISSING_STRING)?;
    }

    f.write_char(DELIMITER)?;
    if let Some(item_rgb) = record.item_rgb() {
        write!(f, "{}", item_rgb)?;
    } else {
        f.write_str(MISSING_STRING)?;
    }

    f.write_char(DELIMITER)?;
    if let (Some(block_sizes), Some(block_starts)) = (record.block_sizes(), record.block_starts()) {
        write!(f, "{}", min(block_sizes.len(), block_starts.len()))?;
    } else {
        f.write_str(MISSING_STRING)?;
    }

    f.write_char(DELIMITER)?;
    if let Some(block_sizes) = record.block_sizes() {
        for i in block_sizes {
            write!(f, "{}", i)?;
            f.write_char(COMMA_SEPARATOR)?;
        }
    } else {
        f.write_str(MISSING_STRING)?;
    }

    f.write_char(DELIMITER)?;
    if let Some(block_starts) = record.block_starts() {
        for i in block_starts {
            write!(f, "{}", i)?;
            f.write_char(COMMA_SEPARATOR)?;
        }
        Ok(())
    } else {
        f.write_str(MISSING_STRING)
    }
}

fn format_optional_fields(
    f: &mut fmt::Formatter<'_>,
    optional_fields: &OptionalFields,
) -> fmt::Result {
    if !optional_fields.is_empty() {
        f.write_char(DELIMITER)?;
        write!(f, "{}", optional_fields)?;
    }

    Ok(())
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
    /// The number of blocks is missing.
    MissingNumBlocks,
    /// The num block is missing.
    InvalidNumBlocks(num::ParseIntError),
    /// The block size is missing.
    MissingBlockSizes,
    /// block size invalid.
    InvalidBlockSize(num::ParseIntError),
    /// the block start.
    MissingBlockStarts,
    /// Failed to parse an element from block start.
    InvalidBlockStarts(num::ParseIntError),
    /// The name is invalid.
    InvalidName(name::ParseError),
    /// The score is missing.
    MissingScore,
    /// The score is invalid.
    InvalidScore(score::ParseError),
    /// The strand is missing.
    MissingStrand,
    /// The strand is invalid.
    InvalidStrand(strand::ParseError),
    InvalidColor(color::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingReferenceSequenceName => f.write_str("missing reference sequence name"),
            Self::MissingStartPosition => f.write_str("missing start position"),
            Self::InvalidStartPosition => f.write_str("invalid start position"),
            Self::MissingEndPosition => f.write_str("missing end position"),
            Self::InvalidEndPosition(e) => write!(f, "invalid end position: {}", e),
            Self::MissingName => f.write_str("missing name"),
            Self::InvalidName(e) => write!(f, "invalid name: {}", e),
            Self::MissingScore => f.write_str("missing score"),
            Self::InvalidScore(e) => write!(f, "invalid score: {}", e),
            Self::MissingStrand => f.write_str("missing strand"),
            Self::InvalidStrand(e) => write!(f, "invalid strand: {}", e),
            Self::MissingNumBlocks => f.write_str("invalid numBlocks"),
            Self::MissingBlockStarts => f.write_str("invalid strand"),
            Self::MissingBlockSizes => f.write_str("missing block size"),
            Self::InvalidBlockSize(e) => write!(f, "invalid block size: {}", e),
            Self::InvalidBlockStarts(e) => write!(f, "invalid block start: {}", e),
            Self::InvalidColor(e) => write!(f, "invalid color: {}", e),
            Self::InvalidNumBlocks(e) => write!(f, "invalid num blocks: {}", e),
        }
    }
}

impl FromStr for Record<3> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        let standard_fields = parse_bed_3_fields(&mut fields)?;
        let optional_fields = parse_optional_fields(&mut fields);
        Ok(Self::new(standard_fields, optional_fields))
    }
}

impl FromStr for Record<4> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        let standard_fields = parse_bed_4_fields(&mut fields)?;
        let optional_fields = parse_optional_fields(&mut fields);
        Ok(Self::new(standard_fields, optional_fields))
    }
}

impl FromStr for Record<5> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        let standard_fields = parse_bed_5_fields(&mut fields)?;
        let optional_fields = parse_optional_fields(&mut fields);
        Ok(Self::new(standard_fields, optional_fields))
    }
}

impl FromStr for Record<6> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        let standard_fields = parse_bed_6_fields(&mut fields)?;
        let optional_fields = parse_optional_fields(&mut fields);
        Ok(Self::new(standard_fields, optional_fields))
    }
}

impl FromStr for Record<12> {
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

fn parse_bed_12_fields<'a, I>(fields: &mut I) -> Result<StandardFields, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    let mut standard_fields = parse_bed_6_fields(fields)?;

    standard_fields.thick_start = Some(
        fields
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
            })?,
    );
    standard_fields.thick_end = Some(
        fields
            .next()
            .ok_or(ParseError::MissingEndPosition)
            .and_then(|s| s.parse().map_err(ParseError::InvalidEndPosition))?,
    );
    standard_fields.item_rgb = parse_color(fields)?;

    let num_blocks: usize = fields
        .next()
        .ok_or(ParseError::MissingNumBlocks)
        .and_then(|s| s.parse().map_err(ParseError::InvalidNumBlocks))?;
    let mut block_sizes: Vec<usize> = Vec::with_capacity(num_blocks);
    let mut block_starts: Vec<Position> = Vec::with_capacity(num_blocks);
    let mut block_size_split = fields
        .next()
        .ok_or(ParseError::MissingBlockSizes)?
        .split(",");
    let mut block_start_split = fields
        .next()
        .ok_or(ParseError::MissingBlockStarts)?
        .split(",");
    while let (Some(size), Some(start)) = (block_size_split.next(), block_start_split.next()) {
        block_sizes.push(size.parse().map_err(ParseError::InvalidBlockSize)?);
        block_starts.push(start.parse().map_err(ParseError::InvalidBlockSize)?);
    }

    standard_fields.block_starts = Some(block_starts);
    standard_fields.block_sizes = Some(block_sizes);

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
        .ok_or(ParseError::MissingName)
        .and_then(|s| match s {
            MISSING_STRING => Ok(None),
            _ => s.parse().map(Some).map_err(ParseError::InvalidName),
        })
}

fn parse_score<'a, I>(fields: &mut I) -> Result<Option<Score>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingName)
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

fn parse_color<'a, I>(fields: &mut I) -> Result<Option<Color>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingStrand)
        .and_then(|s| match s {
            MISSING_STRING => Ok(None),
            _ => s.parse().map(Some).map_err(ParseError::InvalidColor),
        })
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
    fn test_fmt_for_record_3() -> Result<(), noodles_core::position::TryFromIntError> {
        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;

        let standard_fields = StandardFields::new("sq0", start, end);
        let record: Record<3> = Record::new(standard_fields, OptionalFields::default());
        assert_eq!(record.to_string(), "sq0\t7\t13");

        let standard_fields = StandardFields::new("sq0", start, end);
        let record: Record<3> = Record::new(
            standard_fields,
            OptionalFields::from(vec![String::from("ndls")]),
        );
        assert_eq!(record.to_string(), "sq0\t7\t13\tndls");

        Ok(())
    }

    #[test]
    fn test_fmt_for_record_4() -> Result<(), Box<dyn std::error::Error>> {
        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;

        let standard_fields = StandardFields::new("sq0", start, end);
        let record: Record<4> = Record::new(standard_fields, OptionalFields::default());
        assert_eq!(record.to_string(), "sq0\t7\t13\t.");

        let mut standard_fields = StandardFields::new("sq0", start, end);
        standard_fields.name = "ndls1".parse().map(Some)?;
        let record: Record<4> = Record::new(standard_fields, OptionalFields::default());
        assert_eq!(record.to_string(), "sq0\t7\t13\tndls1");

        let standard_fields = StandardFields::new("sq0", start, end);
        let record: Record<4> = Record::new(
            standard_fields,
            OptionalFields::from(vec![String::from("ndls")]),
        );
        assert_eq!(record.to_string(), "sq0\t7\t13\t.\tndls");

        Ok(())
    }

    #[test]
    fn test_fmt_for_record_5() -> Result<(), Box<dyn std::error::Error>> {
        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;

        let standard_fields = StandardFields::new("sq0", start, end);
        let record: Record<5> = Record::new(standard_fields, OptionalFields::default());
        assert_eq!(record.to_string(), "sq0\t7\t13\t.\t0");

        let mut standard_fields = StandardFields::new("sq0", start, end);
        standard_fields.score = Score::try_from(21).map(Some)?;
        let record: Record<5> = Record::new(standard_fields, OptionalFields::default());
        assert_eq!(record.to_string(), "sq0\t7\t13\t.\t21");

        let standard_fields = StandardFields::new("sq0", start, end);
        let record: Record<5> = Record::new(
            standard_fields,
            OptionalFields::from(vec![String::from("ndls")]),
        );
        assert_eq!(record.to_string(), "sq0\t7\t13\t.\t0\tndls");

        Ok(())
    }

    #[test]
    fn test_fmt_for_record_6() -> Result<(), noodles_core::position::TryFromIntError> {
        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;

        let standard_fields = StandardFields::new("sq0", start, end);
        let record: Record<6> = Record::new(standard_fields, OptionalFields::default());
        assert_eq!(record.to_string(), "sq0\t7\t13\t.\t0\t.");

        let mut standard_fields = StandardFields::new("sq0", start, end);
        standard_fields.strand = Some(Strand::Forward);
        let record: Record<6> = Record::new(standard_fields, OptionalFields::default());
        assert_eq!(record.to_string(), "sq0\t7\t13\t.\t0\t+");

        let standard_fields = StandardFields::new("sq0", start, end);
        let record: Record<6> = Record::new(
            standard_fields,
            OptionalFields::from(vec![String::from("ndls")]),
        );
        assert_eq!(record.to_string(), "sq0\t7\t13\t.\t0\t.\tndls");

        Ok(())
    }

    #[test]
    fn test_from_str_for_record_3() -> Result<(), noodles_core::position::TryFromIntError> {
        let actual = "sq0\t7\t13".parse::<Record<3>>();

        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        let standard_fields = StandardFields::new("sq0", start, end);
        let expected = Ok(Record::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_from_str_for_record_4() -> Result<(), Box<dyn std::error::Error>> {
        let actual = "sq0\t7\t13\tndls1".parse::<Record<4>>();

        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        let mut standard_fields = StandardFields::new("sq0", start, end);
        standard_fields.name = "ndls1".parse().map(Some)?;

        let expected = Ok(Record::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_from_str_for_record_5() -> Result<(), Box<dyn std::error::Error>> {
        let actual = "sq0\t7\t13\t.\t21".parse::<Record<5>>();

        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        let mut standard_fields = StandardFields::new("sq0", start, end);
        standard_fields.score = Score::try_from(21).map(Some)?;

        let expected = Ok(Record::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_from_str_for_record_6() -> Result<(), noodles_core::position::TryFromIntError> {
        let actual = "sq0\t7\t13\t.\t0\t+".parse::<Record<6>>();

        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        let mut standard_fields = StandardFields::new("sq0", start, end);
        standard_fields.strand = Some(Strand::Forward);

        let expected = Ok(Record::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);

        Ok(())
    }
}
