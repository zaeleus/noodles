//! BED record and fields.

pub mod score;
pub mod strand;

pub use self::{score::Score, strand::Strand};

use std::{
    error,
    fmt::{self, Write},
    num,
    ops::Deref,
    str::FromStr,
};

const DELIMITER: char = '\t';
const MISSING_STRING: &str = ".";
const MISSING_NUMBER: &str = "0";

#[derive(Clone, Debug, Eq, PartialEq)]
struct StandardFields {
    reference_sequence_name: String,
    start_position: u64,
    end_position: u64,
    name: Option<String>,
    score: Option<Score>,
    strand: Option<Strand>,
}

impl StandardFields {
    fn new<N>(reference_sequence_name: N, start_position: u64, end_position: u64) -> Self
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

impl BedN<4> for Record<4> {}
impl BedN<4> for Record<5> {}
impl BedN<4> for Record<6> {}

impl BedN<5> for Record<5> {}
impl BedN<5> for Record<6> {}

impl BedN<6> for Record<6> {}

impl<const N: u8> Record<N>
where
    Self: BedN<3>,
{
    fn new(standard_fields: StandardFields, optional_fields: OptionalFields) -> Self {
        Self {
            standard_fields,
            optional_fields,
        }
    }

    /// Returns the reference sequence name (`chrom`).
    pub fn reference_sequence_name(&self) -> &str {
        &self.standard_fields.reference_sequence_name
    }

    /// Returns the feature start position (`chromStart`).
    pub fn start_position(&self) -> u64 {
        self.standard_fields.start_position
    }

    /// Returns the feature end position (`chromEnd`).
    pub fn end_position(&self) -> u64 {
        self.standard_fields.end_position
    }

    /// Returns the list of raw optional fields.
    pub fn optional_fields(&self) -> &OptionalFields {
        &self.optional_fields
    }
}

impl<const N: u8> Record<N>
where
    Self: BedN<4>,
{
    /// Returns the feature name (`name`).
    pub fn name(&self) -> Option<&str> {
        self.standard_fields.name.as_deref()
    }
}

impl<const N: u8> Record<N>
where
    Self: BedN<5>,
{
    /// Returns the score (`score`).
    pub fn score(&self) -> Option<Score> {
        self.standard_fields.score
    }
}

impl<const N: u8> Record<N>
where
    Self: BedN<6>,
{
    /// Returns the feature strand (`strand`).
    pub fn strand(&self) -> Option<Strand> {
        self.standard_fields.strand
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

fn format_bed_3_fields<const N: u8>(f: &mut fmt::Formatter<'_>, record: &Record<N>) -> fmt::Result
where
    Record<N>: BedN<3>,
{
    write!(
        f,
        "{}{}{}{}{}",
        record.reference_sequence_name(),
        DELIMITER,
        record.start_position(),
        DELIMITER,
        record.end_position()
    )
}

fn format_bed_4_fields(f: &mut fmt::Formatter<'_>, record: &Record<4>) -> fmt::Result {
    format_bed_3_fields(f, record)?;

    write!(
        f,
        "{}{}",
        DELIMITER,
        record.name().unwrap_or(MISSING_STRING)
    )
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
    InvalidStartPosition(num::ParseIntError),
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
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingReferenceSequenceName => f.write_str("missing reference sequence name"),
            Self::MissingStartPosition => f.write_str("missing start position"),
            Self::InvalidStartPosition(e) => write!(f, "invalid start position: {}", e),
            Self::MissingEndPosition => f.write_str("missing end position"),
            Self::InvalidEndPosition(e) => write!(f, "invalid end position: {}", e),
            Self::MissingName => f.write_str("missing name"),
            Self::MissingScore => f.write_str("missing score"),
            Self::InvalidScore(e) => write!(f, "invalid score: {}", e),
            Self::MissingStrand => f.write_str("missing strand"),
            Self::InvalidStrand(e) => write!(f, "invalid strand: {}", e),
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
        .and_then(|s| s.parse().map_err(ParseError::InvalidStartPosition))?;

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

fn parse_name<'a, I>(fields: &mut I) -> Result<Option<String>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingName)
        .map(|s| match s {
            MISSING_STRING => None,
            _ => Some(s.into()),
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
    fn test_fmt_for_record_3() {
        let standard_fields = StandardFields::new("sq0", 8, 13);
        let record: Record<3> = Record::new(standard_fields, OptionalFields::default());
        assert_eq!(record.to_string(), "sq0\t8\t13");

        let standard_fields = StandardFields::new("sq0", 8, 13);
        let record: Record<3> = Record::new(
            standard_fields,
            OptionalFields::from(vec![String::from("ndls")]),
        );
        assert_eq!(record.to_string(), "sq0\t8\t13\tndls");
    }

    #[test]
    fn test_fmt_for_record_4() {
        let standard_fields = StandardFields::new("sq0", 8, 13);
        let record: Record<4> = Record::new(standard_fields, OptionalFields::default());
        assert_eq!(record.to_string(), "sq0\t8\t13\t.");

        let mut standard_fields = StandardFields::new("sq0", 8, 13);
        standard_fields.name = Some(String::from("ndls1"));
        let record: Record<4> = Record::new(standard_fields, OptionalFields::default());
        assert_eq!(record.to_string(), "sq0\t8\t13\tndls1");

        let standard_fields = StandardFields::new("sq0", 8, 13);
        let record: Record<4> = Record::new(
            standard_fields,
            OptionalFields::from(vec![String::from("ndls")]),
        );
        assert_eq!(record.to_string(), "sq0\t8\t13\t.\tndls");
    }

    #[test]
    fn test_from_str_for_record_3() {
        let actual = "sq0\t8\t13".parse::<Record<3>>();

        let standard_fields = StandardFields::new("sq0", 8, 13);
        let expected = Ok(Record::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_from_str_for_record_4() {
        let actual = "sq0\t8\t13\tndls1".parse::<Record<4>>();

        let mut standard_fields = StandardFields::new("sq0", 8, 13);
        standard_fields.name = Some(String::from("ndls1"));

        let expected = Ok(Record::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_from_str_for_record_5() -> Result<(), score::TryFromIntError> {
        let actual = "sq0\t8\t13\t.\t21".parse::<Record<5>>();

        let mut standard_fields = StandardFields::new("sq0", 8, 13);
        standard_fields.score = Score::try_from(21).map(Some)?;

        let expected = Ok(Record::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_from_str_for_record_6() {
        let actual = "sq0\t8\t13\t.\t0\t+".parse::<Record<6>>();

        let mut standard_fields = StandardFields::new("sq0", 8, 13);
        standard_fields.strand = Some(Strand::Forward);

        let expected = Ok(Record::new(standard_fields, OptionalFields::default()));

        assert_eq!(actual, expected);
    }
}
