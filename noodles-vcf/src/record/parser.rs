use std::{error, fmt};

use super::{
    alternate_bases, chromosome, filters, genotypes, ids, info, position, quality_score,
    reference_bases, Field, Filters, Genotypes, Info, QualityScore, Record, FIELD_DELIMITER,
    MISSING_FIELD,
};
use crate::Header;

/// An error returned when a raw VCF record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is missing.
    MissingField(Field),
    /// The chromosome is invalid.
    InvalidChromosome(chromosome::ParseError),
    /// The position is invalid.
    InvalidPosition(position::ParseError),
    /// The ID is invalid.
    InvalidIds(ids::ParseError),
    /// The reference bases are invalid.
    InvalidReferenceBases(reference_bases::ParseError),
    /// The alternate bases are invalid.
    InvalidAlternateBases(alternate_bases::ParseError),
    /// The quality score is invalid.
    InvalidQualityScore(quality_score::ParseError),
    /// A filter is invalid.
    InvalidFilters(filters::ParseError),
    /// The info is invalid.
    InvalidInfo(info::ParseError),
    /// A genotype is invalid.
    InvalidGenotypes(genotypes::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(field) => write!(f, "missing field: {}", field),
            Self::InvalidChromosome(e) => write!(f, "invalid chromosome: {}", e),
            Self::InvalidPosition(e) => write!(f, "invalid position: {}", e),
            Self::InvalidIds(e) => write!(f, "invalid IDs: {}", e),
            Self::InvalidReferenceBases(e) => write!(f, "invalid reference bases: {}", e),
            Self::InvalidAlternateBases(e) => write!(f, "invalid alternate bases: {}", e),
            Self::InvalidQualityScore(e) => write!(f, "invalid quality score: {}", e),
            Self::InvalidFilters(e) => write!(f, "invalid filters: {}", e),
            Self::InvalidInfo(e) => write!(f, "invalid info: {}", e),
            Self::InvalidGenotypes(e) => write!(f, "invalid genotypes: {}", e),
        }
    }
}

pub fn parse(s: &str, header: &Header) -> Result<Record, ParseError> {
    const MAX_FIELDS: usize = 9;

    let mut fields = s.splitn(MAX_FIELDS, FIELD_DELIMITER);

    let chrom = parse_string(&mut fields, Field::Chromosome)
        .and_then(|s| s.parse().map_err(ParseError::InvalidChromosome))?;

    let pos = parse_string(&mut fields, Field::Position)
        .and_then(|s| s.parse().map_err(ParseError::InvalidPosition))?;

    let ids = parse_string(&mut fields, Field::Ids)
        .and_then(|s| s.parse().map_err(ParseError::InvalidIds))?;

    let r#ref = parse_string(&mut fields, Field::ReferenceBases)
        .and_then(|s| s.parse().map_err(ParseError::InvalidReferenceBases))?;

    let alt = parse_string(&mut fields, Field::AlternateBases)
        .and_then(|s| s.parse().map_err(ParseError::InvalidAlternateBases))?;

    let qual = parse_quality_score(&mut fields)?;
    let filter = parse_filters(&mut fields)?;

    let info = parse_string(&mut fields, Field::Info)
        .and_then(|s| Info::try_from_str(s, header.infos()).map_err(ParseError::InvalidInfo))?;

    let genotypes = if let Some(s) = fields.next() {
        Genotypes::parse(s, header).map_err(ParseError::InvalidGenotypes)?
    } else {
        Genotypes::default()
    };

    Ok(Record {
        chromosome: chrom,
        position: pos,
        ids,
        reference_bases: r#ref,
        alternate_bases: alt,
        quality_score: qual,
        filters: filter,
        info,
        genotypes,
    })
}

fn parse_string<'a, I>(fields: &mut I, field: Field) -> Result<&'a str, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields.next().ok_or(ParseError::MissingField(field))
}

fn parse_quality_score<'a, I>(fields: &mut I) -> Result<Option<QualityScore>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    parse_string(fields, Field::QualityScore).and_then(|s| match s {
        MISSING_FIELD => Ok(None),
        _ => s.parse().map(Some).map_err(ParseError::InvalidQualityScore),
    })
}

fn parse_filters<'a, I>(fields: &mut I) -> Result<Option<Filters>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    parse_string(fields, Field::Filters).and_then(|s| match s {
        MISSING_FIELD => Ok(None),
        _ => s.parse().map(Some).map_err(ParseError::InvalidFilters),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), Box<dyn std::error::Error>> {
        use alternate_bases::Allele;
        use chromosome::Chromosome;
        use reference_bases::Base;

        let s = "chr1\t13\tnd0\tATCG\tA\t5.8\tPASS\tSVTYPE=DEL";
        let record: Record = s.parse()?;

        assert!(matches!(record.chromosome(), Chromosome::Name(name) if name == "chr1"));

        assert_eq!(i32::from(record.position()), 13);

        let ids = record.ids();
        assert_eq!(ids.len(), 1);
        let id: ids::Id = "nd0".parse()?;
        assert!(ids.contains(&id));

        let reference_bases = [Base::A, Base::T, Base::C, Base::G];
        assert_eq!(&record.reference_bases()[..], &reference_bases[..]);

        let alternate_bases = [Allele::Bases(vec![Base::A])];
        assert_eq!(&record.alternate_bases()[..], &alternate_bases[..]);

        assert_eq!(record.quality_score().map(f32::from), Some(5.8));
        assert_eq!(record.filters(), Some(&Filters::Pass));
        assert_eq!(record.info().len(), 1);
        assert!(record.genotypes().is_empty());

        Ok(())
    }

    #[test]
    fn test_from_str_with_genotype_info() -> Result<(), Box<dyn std::error::Error>> {
        use self::genotypes::{
            genotype::{field::Value, Field},
            Genotype,
        };
        use crate::header::format::Key;

        let s = "chr1\t13\tnd0\tATCG\tA\t5.8\tPASS\tSVTYPE=DEL\tGT:GQ\t0|1:13";
        let record: Record = s.parse()?;

        let keys = genotypes::Keys::try_from(vec![Key::Genotype, Key::ConditionalGenotypeQuality])?;
        let genotypes = vec![Genotype::try_from(vec![
            Field::new(Key::Genotype, Some(Value::String(String::from("0|1")))),
            Field::new(Key::ConditionalGenotypeQuality, Some(Value::Integer(13))),
        ])?];

        let actual = record.genotypes();
        let expected = Genotypes::new(keys, genotypes);

        assert_eq!(actual, &expected);

        Ok(())
    }
}
