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

const NULL_FIELD: &str = ".";
const FIELD_DELIMITER: char = '\t';
const MAX_FIELDS: usize = 9;

/// A GFF record.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    reference_sequence_name: String,
    source: String,
    feature: String,
    start: i32,
    end: i32,
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
    /// assert_eq!(record.feature(), ".");
    /// ```
    pub fn feature(&self) -> &str {
        &self.feature
    }

    /// Returns the start position of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let record = gff::Record::default();
    /// assert_eq!(record.start(), 1);
    /// ```
    pub fn start(&self) -> i32 {
        self.start
    }

    /// Returns the end position of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let record = gff::Record::default();
    /// assert_eq!(record.end(), 1);
    /// ```
    pub fn end(&self) -> i32 {
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
    /// assert!(record.attributes().entries().is_empty());
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

/// An error returned when a raw GFF record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A field is missing.
    MissingField(Field),
    /// A field is empty.
    EmptyField(Field),
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
    /// The attributes are invalid.
    InvalidAttributes(attributes::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.splitn(MAX_FIELDS, FIELD_DELIMITER);

        let reference_sequence_name =
            parse_string(&mut fields, Field::ReferenceSequenceName).map(|s| s.into())?;
        let source = parse_string(&mut fields, Field::Source).map(|s| s.into())?;
        let feature = parse_string(&mut fields, Field::Feature).map(|s| s.into())?;

        let start = parse_string(&mut fields, Field::Start)
            .and_then(|s| s.parse().map_err(ParseError::InvalidStart))?;

        let end = parse_string(&mut fields, Field::End)
            .and_then(|s| s.parse().map_err(ParseError::InvalidEnd))?;

        let score = parse_string(&mut fields, Field::Score).and_then(|s| {
            if s == NULL_FIELD {
                Ok(None)
            } else {
                s.parse().map(Some).map_err(ParseError::InvalidScore)
            }
        })?;

        let strand = parse_string(&mut fields, Field::Strand)
            .and_then(|s| s.parse().map_err(ParseError::InvalidStrand))?;

        let phase = parse_string(&mut fields, Field::Phase).and_then(|s| {
            if s == NULL_FIELD {
                Ok(None)
            } else {
                s.parse().map(Some).map_err(ParseError::InvalidPhase)
            }
        })?;

        let attributes = match fields.next() {
            Some(s) => s.parse().map_err(ParseError::InvalidAttributes)?,
            None => Attributes::default(),
        };

        Ok(Record {
            reference_sequence_name,
            source,
            feature,
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
    fields.next().ok_or_else(|| ParseError::MissingField(field))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let s = "sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id=ndls0;gene_name=gene0";
        let record = s.parse::<Record>()?;

        assert_eq!(record.reference_sequence_name(), "sq0");
        assert_eq!(record.source(), "NOODLES");
        assert_eq!(record.feature(), "gene");
        assert_eq!(record.start(), 8);
        assert_eq!(record.end(), 13);
        assert_eq!(record.score(), None);
        assert_eq!(record.strand(), Strand::Forward);
        assert_eq!(record.phase(), None);

        assert_eq!(
            record.attributes(),
            &Attributes::from(vec![
                attributes::Entry::new(String::from("gene_id"), String::from("ndls0")),
                attributes::Entry::new(String::from("gene_name"), String::from("gene0")),
            ])
        );

        Ok(())
    }
}
