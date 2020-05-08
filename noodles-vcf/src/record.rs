mod chromosome;
mod field;
mod filter_status;
mod reference_bases;

pub use self::{
    chromosome::Chromosome, field::Field, filter_status::FilterStatus,
    reference_bases::ReferenceBases,
};

use std::{error, fmt, str::FromStr};

pub(crate) const MISSING_FIELD: &str = ".";
const FIELD_DELIMITER: char = '\t';

#[derive(Debug)]
pub struct Record {
    chromosome: Chromosome,
    position: i32,
    id: String,
    reference_bases: ReferenceBases,
    alternate_bases: String,
    quality_score: f32,
    filter_status: FilterStatus,
    information: String,
}

impl Record {
    pub fn chromosome(&self) -> &Chromosome {
        &self.chromosome
    }

    pub fn position(&self) -> i32 {
        self.position
    }

    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn reference_bases(&self) -> &ReferenceBases {
        &self.reference_bases
    }

    pub fn alternate_bases(&self) -> &str {
        &self.alternate_bases
    }

    pub fn quality_score(&self) -> f32 {
        self.quality_score
    }

    pub fn filter_status(&self) -> &FilterStatus {
        &self.filter_status
    }

    pub fn information(&self) -> &str {
        &self.information
    }
}

#[derive(Debug)]
pub enum ParseError {
    Missing(Field),
    Invalid(Field, Box<dyn error::Error>),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Missing(field) => write!(f, "missing field: {}", field),
            Self::Invalid(field, message) => write!(f, "invalid {} field: {}", field, message),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(FIELD_DELIMITER);

        let chrom = parse_string(&mut fields, Field::Chromosome).and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::Invalid(Field::Chromosome, Box::new(e)))
        })?;

        let pos = parse_i32(&mut fields, Field::Position)?;
        let id = parse_string(&mut fields, Field::Id)?;

        let r#ref = parse_string(&mut fields, Field::ReferenceBases).and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::Invalid(Field::ReferenceBases, Box::new(e)))
        })?;

        let alt = parse_string(&mut fields, Field::AlternateBases)?;
        let qual = parse_f32(&mut fields, Field::QualityScore)?;

        let filter = parse_string(&mut fields, Field::FilterStatus).and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::Invalid(Field::FilterStatus, Box::new(e)))
        })?;

        let info = parse_string(&mut fields, Field::Information)?;

        Ok(Self {
            chromosome: chrom,
            position: pos,
            id: id.into(),
            reference_bases: r#ref,
            alternate_bases: alt.into(),
            quality_score: qual,
            filter_status: filter,
            information: info.into(),
        })
    }
}

fn parse_string<'a, I>(fields: &mut I, field: Field) -> Result<&'a str, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields.next().ok_or_else(|| ParseError::Missing(field))
}

fn parse_i32<'a, I>(fields: &mut I, field: Field) -> Result<i32, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    parse_string(fields, field).and_then(|s| {
        s.parse()
            .map_err(|e| ParseError::Invalid(field, Box::new(e)))
    })
}

fn parse_f32<'a, I>(fields: &mut I, field: Field) -> Result<f32, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    parse_string(fields, field).and_then(|s| {
        s.parse()
            .map_err(|e| ParseError::Invalid(field, Box::new(e)))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        use reference_bases::Base;

        let s = "chr1\t13\tr0\tATCG\tA\t5.8\tPASS\tSVTYPE=DEL";
        let record: Record = s.parse()?;

        assert!(matches!(record.chromosome(), Chromosome::Name(name) if name == "chr1"));

        assert_eq!(record.position(), 13);
        assert_eq!(record.id(), "r0");

        let reference_bases = [Base::A, Base::T, Base::C, Base::G];
        assert_eq!(&record.reference_bases()[..], &reference_bases[..]);

        assert_eq!(record.alternate_bases(), "A");
        assert_eq!(record.quality_score(), 5.8);
        assert_eq!(record.filter_status(), &FilterStatus::Pass);
        assert_eq!(record.information(), "SVTYPE=DEL");

        Ok(())
    }
}
