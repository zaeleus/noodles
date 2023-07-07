//! GFF record and fields.

pub mod attributes;
mod builder;
mod field;
mod phase;
mod strand;

pub use self::{
    attributes::Attributes, builder::Builder, field::Field, phase::Phase, strand::Strand,
};

use std::{error, fmt, num, str::FromStr};

use noodles_core::Position;

pub(crate) const MISSING_FIELD: &str = ".";
const FIELD_DELIMITER: char = '\t';
const MAX_FIELDS: usize = 9;

/// A GFF record.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    reference_sequence_name: String,
    source: String,
    ty: String,
    start: Position,
    end: Position,
    score: Option<f32>,
    strand: Strand,
    phase: Option<Phase>,
    attributes: Attributes,
}

impl Record {
    /// Returns a builder to create a record from each of its fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    ///
    /// let record = gff::Record::builder()
    ///     .set_reference_sequence_name(String::from("sq0"))
    ///     .build();
    ///
    /// assert_eq!(record.reference_sequence_name(), "sq0");
    /// ```
    pub fn builder() -> Builder {
        Builder::new()
    }

    /// Returns the reference sequence name of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let record = gff::Record::default();
    /// assert_eq!(record.reference_sequence_name(), ".");
    /// ```
    pub fn reference_sequence_name(&self) -> &str {
        &self.reference_sequence_name
    }

    /// Returns the source of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let record = gff::Record::default();
    /// assert_eq!(record.source(), ".");
    /// ```
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Returns the feature type of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let record = gff::Record::default();
    /// assert_eq!(record.ty(), ".");
    /// ```
    pub fn ty(&self) -> &str {
        &self.ty
    }

    /// Returns the start position of the record.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_gff as gff;
    /// let record = gff::Record::default();
    /// assert_eq!(record.start(), Position::MIN);
    /// ```
    pub fn start(&self) -> Position {
        self.start
    }

    /// Returns the end position of the record.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_gff as gff;
    /// let record = gff::Record::default();
    /// assert_eq!(record.end(), Position::MIN);
    /// ```
    pub fn end(&self) -> Position {
        self.end
    }

    /// Returns the score of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let record = gff::Record::default();
    /// assert!(record.score().is_none());
    /// ```
    pub fn score(&self) -> Option<f32> {
        self.score
    }

    /// Returns the strand of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::{self as gff, record::Strand};
    /// let record = gff::Record::default();
    /// assert_eq!(record.strand(), Strand::None);
    /// ```
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Returns the phase of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let record = gff::Record::default();
    /// assert!(record.phase().is_none());
    /// ```
    pub fn phase(&self) -> Option<Phase> {
        self.phase
    }

    /// Returns the attributes of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let record = gff::Record::default();
    /// assert!(record.attributes().is_empty());
    /// ```
    pub fn attributes(&self) -> &Attributes {
        &self.attributes
    }
}

impl Default for Record {
    fn default() -> Self {
        Builder::new().build()
    }
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{seqid}\t{source}\t{ty}\t{start}\t{end}",
            seqid = self.reference_sequence_name(),
            source = self.source(),
            ty = self.ty(),
            start = self.start(),
            end = self.end(),
        )?;

        if let Some(score) = self.score() {
            write!(f, "\t{score}")?;
        } else {
            write!(f, "\t{MISSING_FIELD}")?;
        }

        write!(f, "\t{}", self.strand())?;

        if let Some(phase) = self.phase() {
            write!(f, "\t{phase}")?;
        } else {
            write!(f, "\t{MISSING_FIELD}")?;
        }

        if self.attributes().is_empty() {
            write!(f, "\t{MISSING_FIELD}")?;
        } else {
            write!(f, "\t{}", self.attributes())?;
        }

        Ok(())
    }
}

/// An error returned when a raw GFF record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A field is missing.
    MissingField(Field),
    /// A field is empty.
    EmptyField(Field),
    /// The reference sequence name is invalid.
    InvalidReferenceSequenceName,
    /// The start is invalid.
    InvalidStart(num::ParseIntError),
    /// The end is invalid.
    InvalidEnd(num::ParseIntError),
    /// The score is invalid.
    InvalidScore(num::ParseFloatError),
    /// The strand is invalid.
    InvalidStrand(strand::ParseError),
    /// The phase is invalid.
    InvalidPhase(phase::ParseError),
    /// The phase is missing.
    ///
    /// The phase is required for CDS features.
    MissingPhase,
    /// The attributes are invalid.
    InvalidAttributes(attributes::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidStart(e) | Self::InvalidEnd(e) => Some(e),
            Self::InvalidScore(e) => Some(e),
            Self::InvalidStrand(e) => Some(e),
            Self::InvalidPhase(e) => Some(e),
            Self::InvalidAttributes(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::MissingField(field) => write!(f, "missing field: {field:?}"),
            Self::EmptyField(field) => write!(f, "empty field: {field:?}"),
            Self::InvalidReferenceSequenceName => write!(f, "invalid reference sequence name"),
            Self::InvalidStart(_) => f.write_str("invalid start"),
            Self::InvalidEnd(_) => f.write_str("invalid end"),
            Self::InvalidScore(_) => f.write_str("invalid score"),
            Self::InvalidStrand(_) => f.write_str("invalid strand"),
            Self::InvalidPhase(_) => f.write_str("invalid phase"),
            Self::MissingPhase => write!(f, "missing phase"),
            Self::InvalidAttributes(_) => f.write_str("invalid attributes"),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.splitn(MAX_FIELDS, FIELD_DELIMITER);

        let reference_sequence_name = parse_string(&mut fields, Field::ReferenceSequenceName)
            .and_then(parse_reference_sequence_name)?;

        let source = parse_string(&mut fields, Field::Source).map(|s| s.into())?;
        let ty = parse_string(&mut fields, Field::Type).map(|s| s.into())?;

        let start = parse_string(&mut fields, Field::Start)
            .and_then(|s| s.parse().map_err(ParseError::InvalidStart))?;

        let end = parse_string(&mut fields, Field::End)
            .and_then(|s| s.parse().map_err(ParseError::InvalidEnd))?;

        let score = parse_string(&mut fields, Field::Score).and_then(|s| {
            if s == MISSING_FIELD {
                Ok(None)
            } else {
                s.parse().map(Some).map_err(ParseError::InvalidScore)
            }
        })?;

        let strand = parse_string(&mut fields, Field::Strand)
            .and_then(|s| s.parse().map_err(ParseError::InvalidStrand))?;

        let phase = parse_string(&mut fields, Field::Phase).and_then(|s| {
            if s == MISSING_FIELD {
                if ty == "CDS" {
                    Err(ParseError::MissingPhase)
                } else {
                    Ok(None)
                }
            } else {
                s.parse().map(Some).map_err(ParseError::InvalidPhase)
            }
        })?;

        let attributes = match fields.next() {
            Some(s) => s.parse().map_err(ParseError::InvalidAttributes)?,
            None => Attributes::default(),
        };

        Ok(Self {
            reference_sequence_name,
            source,
            ty,
            start,
            end,
            score,
            strand,
            phase,
            attributes,
        })
    }
}

fn parse_string<'a, I>(fields: &mut I, field: Field) -> Result<&'a str, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields.next().ok_or(ParseError::MissingField(field))
}

fn parse_reference_sequence_name(s: &str) -> Result<String, ParseError> {
    if s.starts_with('>') {
        Err(ParseError::InvalidReferenceSequenceName)
    } else {
        Ok(s.into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let record = Record::default();
        assert_eq!(record.to_string(), ".\t.\t.\t1\t1\t.\t.\t.\t.");
    }

    #[test]
    fn test_from_str() -> Result<(), Box<dyn std::error::Error>> {
        use self::attributes::field::{Key, Value};

        let s = "sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id=ndls0;gene_name=gene0";
        let record = s.parse::<Record>()?;

        assert_eq!(record.reference_sequence_name(), "sq0");
        assert_eq!(record.source(), "NOODLES");
        assert_eq!(record.ty(), "gene");
        assert_eq!(record.start(), Position::try_from(8)?);
        assert_eq!(record.end(), Position::try_from(13)?);
        assert_eq!(record.score(), None);
        assert_eq!(record.strand(), Strand::Forward);
        assert_eq!(record.phase(), None);

        assert_eq!(
            record.attributes(),
            &[
                (Key::from("gene_id"), Value::from("ndls0")),
                (Key::from("gene_name"), Value::from("gene0")),
            ]
            .into_iter()
            .collect()
        );

        Ok(())
    }

    #[test]
    fn test_from_str_with_cds_feature_and_no_phase() {
        let s = "sq0\tNOODLES\tCDS\t8\t13\t.\t+\t.\tgene_id=ndls0;gene_name=gene0";
        assert_eq!(s.parse::<Record>(), Err(ParseError::MissingPhase));
    }

    #[test]
    fn test_parse_reference_sequence_name() {
        assert_eq!(
            parse_reference_sequence_name("sq0"),
            Ok(String::from("sq0"))
        );

        assert_eq!(
            parse_reference_sequence_name(">sq0"),
            Err(ParseError::InvalidReferenceSequenceName)
        );
    }
}
