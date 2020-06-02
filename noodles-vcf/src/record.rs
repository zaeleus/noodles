pub mod alternate_bases;
mod builder;
mod chromosome;
mod field;
mod filter_status;
pub mod format;
pub mod genotype;
mod id;
pub mod info;
mod quality_score;
mod reference_bases;

pub use self::{
    alternate_bases::AlternateBases, builder::Builder, chromosome::Chromosome, field::Field,
    filter_status::FilterStatus, format::Format, genotype::Genotype, id::Id, info::Info,
    quality_score::QualityScore, reference_bases::ReferenceBases,
};

use std::{error, fmt, num, str::FromStr};

pub(crate) const MISSING_FIELD: &str = ".";
const FIELD_DELIMITER: char = '\t';

#[derive(Debug)]
pub struct Record {
    chromosome: Chromosome,
    position: i32,
    id: Id,
    reference_bases: ReferenceBases,
    alternate_bases: AlternateBases,
    quality_score: QualityScore,
    filter_status: FilterStatus,
    info: Info,
    format: Option<Format>,
    genotypes: Vec<Genotype>,
}

impl Record {
    pub fn builder() -> Builder {
        Builder::new()
    }

    pub fn chromosome(&self) -> &Chromosome {
        &self.chromosome
    }

    pub fn position(&self) -> i32 {
        self.position
    }

    pub fn id(&self) -> &Id {
        &self.id
    }

    pub fn reference_bases(&self) -> &ReferenceBases {
        &self.reference_bases
    }

    pub fn alternate_bases(&self) -> &AlternateBases {
        &self.alternate_bases
    }

    pub fn quality_score(&self) -> &QualityScore {
        &self.quality_score
    }

    pub fn filter_status(&self) -> &FilterStatus {
        &self.filter_status
    }

    pub fn info(&self) -> &Info {
        &self.info
    }

    pub fn format(&self) -> Option<&Format> {
        self.format.as_ref()
    }

    pub fn genotypes(&self) -> &[Genotype] {
        &self.genotypes
    }
}

#[derive(Debug)]
pub enum ParseError {
    MissingField(Field),
    InvalidChromosome(chromosome::ParseError),
    InvalidPosition(num::ParseIntError),
    InvalidId(id::ParseError),
    InvalidReferenceBases(reference_bases::ParseError),
    InvalidAlternateBases(alternate_bases::ParseError),
    InvalidQualityScore(quality_score::ParseError),
    InvalidFilterStatus(filter_status::ParseError),
    InvalidInfo(info::ParseError),
    InvalidFormat(format::ParseError),
    InvalidGenotype(genotype::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(field) => write!(f, "missing field: {}", field),
            Self::InvalidChromosome(e) => write!(f, "{}", e),
            Self::InvalidPosition(e) => write!(f, "{}", e),
            Self::InvalidId(e) => write!(f, "{}", e),
            Self::InvalidReferenceBases(e) => write!(f, "{}", e),
            Self::InvalidAlternateBases(e) => write!(f, "{}", e),
            Self::InvalidQualityScore(e) => write!(f, "{}", e),
            Self::InvalidFilterStatus(e) => write!(f, "{}", e),
            Self::InvalidInfo(e) => write!(f, "{}", e),
            Self::InvalidFormat(e) => write!(f, "{}", e),
            Self::InvalidGenotype(e) => write!(f, "{}", e),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(FIELD_DELIMITER);

        let chrom = parse_string(&mut fields, Field::Chromosome)
            .and_then(|s| s.parse().map_err(ParseError::InvalidChromosome))?;

        let pos = parse_string(&mut fields, Field::Position)
            .and_then(|s| s.parse().map_err(ParseError::InvalidPosition))?;

        let id = parse_string(&mut fields, Field::Id)
            .and_then(|s| s.parse().map_err(ParseError::InvalidId))?;

        let r#ref = parse_string(&mut fields, Field::ReferenceBases)
            .and_then(|s| s.parse().map_err(ParseError::InvalidReferenceBases))?;

        let alt = parse_string(&mut fields, Field::AlternateBases)
            .and_then(|s| s.parse().map_err(ParseError::InvalidAlternateBases))?;

        let qual = parse_string(&mut fields, Field::QualityScore)
            .and_then(|s| s.parse().map_err(ParseError::InvalidQualityScore))?;

        let filter = parse_string(&mut fields, Field::FilterStatus)
            .and_then(|s| s.parse().map_err(ParseError::InvalidFilterStatus))?;

        let info = parse_string(&mut fields, Field::Info)
            .and_then(|s| s.parse().map_err(ParseError::InvalidInfo))?;

        let format = match fields.next() {
            Some(s) => s.parse().map(Some).map_err(ParseError::InvalidFormat)?,
            None => None,
        };

        let genotypes = format
            .as_ref()
            .map(|f| {
                fields
                    .map(|s| Genotype::from_str_format(s, f))
                    .collect::<Result<_, _>>()
                    .map_err(ParseError::InvalidGenotype)
            })
            .unwrap_or_else(|| Ok(Vec::new()))?;

        Ok(Self {
            chromosome: chrom,
            position: pos,
            id,
            reference_bases: r#ref,
            alternate_bases: alt,
            quality_score: qual,
            filter_status: filter,
            info,
            format,
            genotypes,
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
        use alternate_bases::Allele;
        use reference_bases::Base;

        let s = "chr1\t13\tr0\tATCG\tA\t5.8\tPASS\tSVTYPE=DEL";
        let record: Record = s.parse()?;

        assert!(matches!(record.chromosome(), Chromosome::Name(name) if name == "chr1"));

        assert_eq!(record.position(), 13);
        assert!(record.id().is_some());

        let reference_bases = [Base::A, Base::T, Base::C, Base::G];
        assert_eq!(&record.reference_bases()[..], &reference_bases[..]);

        let alternate_bases = [Allele::Bases(vec![Base::A])];
        assert_eq!(&record.alternate_bases()[..], &alternate_bases[..]);

        assert_eq!(**record.quality_score(), Some(5.8));
        assert_eq!(record.filter_status(), &FilterStatus::Pass);
        assert_eq!(record.info().len(), 1);
        assert!(record.format().is_none());
        assert!(record.genotypes().is_empty());

        Ok(())
    }

    #[test]
    fn test_from_str_with_genotype_info() -> Result<(), ParseError> {
        let s = "chr1\t13\tr0\tATCG\tA\t5.8\tPASS\tSVTYPE=DEL\tGT:GQ\t0|1:13";
        let record: Record = s.parse()?;

        let expected = [
            genotype::field::Key::Genotype,
            genotype::field::Key::ConditionalGenotypeQuality,
        ];

        assert_eq!(record.format().map(|f| &f[..]), Some(&expected[..]));

        let genotypes = record.genotypes();
        let expected = vec![
            genotype::Field::new(
                genotype::field::Key::Genotype,
                genotype::field::Value::String(String::from("0|1")),
            ),
            genotype::Field::new(
                genotype::field::Key::ConditionalGenotypeQuality,
                genotype::field::Value::Integer(13),
            ),
        ];

        assert_eq!(genotypes.len(), 1);
        assert_eq!(&genotypes[0][..], &expected[..]);

        Ok(())
    }
}
