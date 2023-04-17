//! VCF header format key.

mod v4_3;
mod v4_4;

use crate::{
    header::{record::value::map::format::Type, FileFormat, Number},
    record::genotypes::keys::Key,
};

pub(crate) fn definition(
    file_format: FileFormat,
    key: &Key,
) -> Option<(Number, Type, &'static str)> {
    match key {
        Key::Standard(k) => match (file_format.major(), file_format.minor()) {
            (4, 4) => v4_4::definition(*k),
            (4, 3) => v4_3::definition(*k),
            _ => None,
        },
        Key::Other(_) => None,
    }
}

pub(crate) fn number(key: &Key) -> Option<Number> {
    match key {
        Key::Standard(k) => v4_3::definition(*k).map(|(number, _, _)| number),
        Key::Other(_) => None,
    }
}

pub(crate) fn ty(key: &Key) -> Option<Type> {
    match key {
        Key::Standard(k) => v4_3::definition(*k).map(|(_, ty, _)| ty),
        Key::Other(_) => None,
    }
}

pub(crate) fn description(key: &Key) -> Option<&str> {
    match key {
        Key::Standard(k) => v4_3::definition(*k).map(|(_, _, description)| description),
        Key::Other(_) => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::genotypes::keys::key;

    #[test]
    fn test_number() -> Result<(), key::ParseError> {
        assert_eq!(number(&key::READ_DEPTHS), Some(Number::R));
        assert_eq!(number(&key::FORWARD_STRAND_READ_DEPTHS), Some(Number::R));
        assert_eq!(number(&key::REVERSE_STRAND_READ_DEPTHS), Some(Number::R));
        assert_eq!(number(&key::READ_DEPTH), Some(Number::Count(1)));
        assert_eq!(
            number(&key::EXPECTED_ALTERNATE_ALLELE_COUNTS),
            Some(Number::A)
        );
        assert_eq!(number(&key::FILTER), Some(Number::Count(1)));
        assert_eq!(number(&key::GENOTYPE_LIKELIHOODS), Some(Number::G));
        assert_eq!(
            number(&key::GENOTYPE_POSTERIOR_PROBABILITIES),
            Some(Number::G)
        );
        assert_eq!(
            number(&key::CONDITIONAL_GENOTYPE_QUALITY),
            Some(Number::Count(1))
        );
        assert_eq!(number(&key::GENOTYPE), Some(Number::Count(1)));
        assert_eq!(number(&key::HAPLOTYPE_QUALITY), Some(Number::Count(2)));
        assert_eq!(number(&key::MAPPING_QUALITY), Some(Number::Count(1)));
        assert_eq!(number(&key::ROUNDED_GENOTYPE_LIKELIHOODS), Some(Number::G));
        assert_eq!(
            number(&key::ROUNDED_GENOTYPE_POSTERIOR_PROBABILITIES),
            Some(Number::G)
        );
        assert_eq!(number(&key::PHASING_QUALITY), Some(Number::Count(1)));
        assert_eq!(number(&key::PHASE_SET), Some(Number::Count(1)));

        assert_eq!(number(&key::GENOTYPE_COPY_NUMBER), Some(Number::Count(1)));
        assert_eq!(
            number(&key::GENOTYPE_COPY_NUMBER_QUALITY),
            Some(Number::Count(1))
        );
        assert_eq!(
            number(&key::GENOTYPE_COPY_NUMBER_LIKELIHOODS),
            Some(Number::G)
        );
        assert_eq!(
            number(&key::GENOTYPE_COPY_NUMBER_POSTERIOR_PROBABILITIES),
            Some(Number::G)
        );
        assert_eq!(
            number(&key::NOVEL_VARIANT_QUALITY_SCORE),
            Some(Number::Count(1))
        );
        assert_eq!(number(&key::HAPLOTYPE_ID), Some(Number::Count(1)));
        assert_eq!(number(&key::ANCESTRAL_HAPLOTYPE_ID), Some(Number::Count(1)));

        assert!(number(&"NDLS".parse()?).is_none());

        Ok(())
    }

    #[test]
    fn test_ty() -> Result<(), key::ParseError> {
        assert_eq!(ty(&key::READ_DEPTHS), Some(Type::Integer));
        assert_eq!(ty(&key::FORWARD_STRAND_READ_DEPTHS), Some(Type::Integer));
        assert_eq!(ty(&key::REVERSE_STRAND_READ_DEPTHS), Some(Type::Integer));
        assert_eq!(ty(&key::READ_DEPTH), Some(Type::Integer));
        assert_eq!(
            ty(&key::EXPECTED_ALTERNATE_ALLELE_COUNTS),
            Some(Type::Integer)
        );
        assert_eq!(ty(&key::FILTER), Some(Type::String));
        assert_eq!(ty(&key::GENOTYPE_LIKELIHOODS), Some(Type::Float));
        assert_eq!(
            ty(&key::GENOTYPE_POSTERIOR_PROBABILITIES),
            Some(Type::Float)
        );
        assert_eq!(ty(&key::CONDITIONAL_GENOTYPE_QUALITY), Some(Type::Integer));
        assert_eq!(ty(&key::GENOTYPE), Some(Type::String));
        assert_eq!(ty(&key::HAPLOTYPE_QUALITY), Some(Type::Integer));
        assert_eq!(ty(&key::MAPPING_QUALITY), Some(Type::Integer));
        assert_eq!(ty(&key::ROUNDED_GENOTYPE_LIKELIHOODS), Some(Type::Integer));
        assert_eq!(
            ty(&key::ROUNDED_GENOTYPE_POSTERIOR_PROBABILITIES),
            Some(Type::Integer)
        );
        assert_eq!(ty(&key::PHASING_QUALITY), Some(Type::Integer));
        assert_eq!(ty(&key::PHASE_SET), Some(Type::Integer));

        assert_eq!(ty(&key::GENOTYPE_COPY_NUMBER), Some(Type::Integer));
        assert_eq!(ty(&key::GENOTYPE_COPY_NUMBER_QUALITY), Some(Type::Float));
        assert_eq!(
            ty(&key::GENOTYPE_COPY_NUMBER_LIKELIHOODS),
            Some(Type::Float)
        );
        assert_eq!(
            ty(&key::GENOTYPE_COPY_NUMBER_POSTERIOR_PROBABILITIES),
            Some(Type::Float)
        );
        assert_eq!(ty(&key::NOVEL_VARIANT_QUALITY_SCORE), Some(Type::Integer));
        assert_eq!(ty(&key::HAPLOTYPE_ID), Some(Type::Integer));
        assert_eq!(ty(&key::ANCESTRAL_HAPLOTYPE_ID), Some(Type::Integer));

        assert!(ty(&"NDLS".parse()?).is_none());

        Ok(())
    }

    #[test]
    fn test_description() -> Result<(), key::ParseError> {
        assert_eq!(
            description(&key::READ_DEPTHS),
            Some("Read depth for each allele")
        );
        assert_eq!(
            description(&key::FORWARD_STRAND_READ_DEPTHS),
            Some("Read depth for each allele on the forward strand")
        );
        assert_eq!(
            description(&key::REVERSE_STRAND_READ_DEPTHS),
            Some("Read depth for each allele on the reverse strand")
        );
        assert_eq!(description(&key::READ_DEPTH), Some("Read depth"));
        assert_eq!(
            description(&key::EXPECTED_ALTERNATE_ALLELE_COUNTS),
            Some("Expected alternate allele counts")
        );
        assert_eq!(
            description(&key::FILTER),
            Some(r#"Filter indicating if this genotype was "called""#)
        );
        assert_eq!(
            description(&key::GENOTYPE_LIKELIHOODS),
            Some("Genotype likelihoods")
        );
        assert_eq!(
            description(&key::GENOTYPE_POSTERIOR_PROBABILITIES),
            Some("Genotype posterior probabilities")
        );
        assert_eq!(
            description(&key::CONDITIONAL_GENOTYPE_QUALITY),
            Some("Conditional genotype quality")
        );
        assert_eq!(description(&key::GENOTYPE), Some("Genotype"));
        assert_eq!(
            description(&key::HAPLOTYPE_QUALITY),
            Some("Haplotype quality")
        );
        assert_eq!(
            description(&key::MAPPING_QUALITY),
            Some("RMS mapping quality")
        );
        assert_eq!(
            description(&key::ROUNDED_GENOTYPE_LIKELIHOODS),
            Some("Phred-scaled genotype likelihoods rounded to the closest integer")
        );
        assert_eq!(
            description(&key::ROUNDED_GENOTYPE_POSTERIOR_PROBABILITIES),
            Some("Phred-scaled genotype posterior probabilities rounded to the closest integer")
        );
        assert_eq!(description(&key::PHASING_QUALITY), Some("Phasing quality"));
        assert_eq!(description(&key::PHASE_SET), Some("Phase set"));

        assert_eq!(
            description(&key::GENOTYPE_COPY_NUMBER),
            Some("Copy number genotype for imprecise events")
        );
        assert_eq!(
            description(&key::GENOTYPE_COPY_NUMBER_QUALITY),
            Some("Copy number genotype quality for imprecise events")
        );
        assert_eq!(
            description(&key::GENOTYPE_COPY_NUMBER_LIKELIHOODS),
            Some("Copy number genotype likelihood for imprecise events")
        );
        assert_eq!(
            description(&key::GENOTYPE_COPY_NUMBER_POSTERIOR_PROBABILITIES),
            Some("Copy number posterior probabilities")
        );
        assert_eq!(
            description(&key::NOVEL_VARIANT_QUALITY_SCORE),
            Some("Phred style probability score that the variant is novel")
        );
        assert_eq!(
            description(&key::HAPLOTYPE_ID),
            Some("Unique haplotype identifier")
        );
        assert_eq!(
            description(&key::ANCESTRAL_HAPLOTYPE_ID),
            Some("Unique identifier of ancestral haplotype")
        );

        assert!(description(&"NDLS".parse()?).is_none());

        Ok(())
    }
}
