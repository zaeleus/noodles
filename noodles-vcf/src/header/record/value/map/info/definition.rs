//! VCF header info key.

mod v4_3;
mod v4_4;

use crate::{
    header::{record::value::map::info::Type, FileFormat, Number},
    record::info::field::Key,
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
    use crate::record::info::field::key;

    #[test]
    fn test_number() -> Result<(), key::ParseError> {
        assert_eq!(number(&key::ANCESTRAL_ALLELE), Some(Number::Count(1)));
        assert_eq!(number(&key::ALLELE_COUNT), Some(Number::A));
        assert_eq!(number(&key::TOTAL_READ_DEPTHS), Some(Number::R));
        assert_eq!(number(&key::FORWARD_STRAND_READ_DEPTHS), Some(Number::R));
        assert_eq!(number(&key::REVERSE_STRAND_READ_DEPTHS), Some(Number::R));
        assert_eq!(number(&key::ALLELE_FREQUENCIES), Some(Number::A));
        assert_eq!(number(&key::TOTAL_ALLELE_COUNT), Some(Number::Count(1)));
        assert_eq!(number(&key::BASE_QUALITY), Some(Number::Count(1)));
        assert_eq!(number(&key::CIGAR), Some(Number::A));
        assert_eq!(number(&key::IS_IN_DB_SNP), Some(Number::Count(0)));
        assert_eq!(number(&key::TOTAL_DEPTH), Some(Number::Count(1)));
        // assert_eq!(number(&END_POSITION), Some(Number::Count(1)));
        assert_eq!(number(&key::IS_IN_HAP_MAP_2), Some(Number::Count(0)));
        assert_eq!(number(&key::IS_IN_HAP_MAP_3), Some(Number::Count(0)));
        assert_eq!(number(&key::MAPPING_QUALITY), Some(Number::Count(1)));
        assert_eq!(
            number(&key::ZERO_MAPPING_QUALITY_COUNT),
            Some(Number::Count(1))
        );
        assert_eq!(
            number(&key::SAMPLES_WITH_DATA_COUNT),
            Some(Number::Count(1))
        );
        assert_eq!(number(&key::STRAND_BIAS), Some(Number::Count(4)));
        assert_eq!(number(&key::IS_SOMATIC_MUTATION), Some(Number::Count(0)));
        assert_eq!(number(&key::IS_VALIDATED), Some(Number::Count(0)));
        assert_eq!(number(&key::IS_IN_1000_GENOMES), Some(Number::Count(0)));

        assert_eq!(number(&key::IS_IMPRECISE), Some(Number::Count(0)));
        assert_eq!(number(&key::IS_NOVEL), Some(Number::Count(0)));
        assert_eq!(number(&key::END_POSITION), Some(Number::Count(1)));
        assert_eq!(number(&key::SV_TYPE), Some(Number::Count(1)));
        assert_eq!(number(&key::SV_LENGTHS), Some(Number::Unknown));
        assert_eq!(
            number(&key::POSITION_CONFIDENCE_INTERVALS),
            Some(Number::Count(2))
        );
        assert_eq!(
            number(&key::END_CONFIDENCE_INTERVALS),
            Some(Number::Count(2))
        );
        assert_eq!(number(&key::MICROHOMOLOGY_LENGTHS), Some(Number::Unknown));
        assert_eq!(number(&key::MICROHOMOLOGY_SEQUENCES), Some(Number::Unknown));
        assert_eq!(number(&key::BREAKPOINT_IDS), Some(Number::Unknown));
        assert_eq!(number(&key::MOBILE_ELEMENT_INFO), Some(Number::Count(4)));
        assert_eq!(
            number(&key::MOBILE_ELEMENT_TRANSDUCTION_INFO),
            Some(Number::Count(4))
        );
        assert_eq!(number(&key::DBV_ID), Some(Number::Count(1)));
        assert_eq!(number(&key::DB_VAR_ID), Some(Number::Count(1)));
        assert_eq!(number(&key::DB_RIP_ID), Some(Number::Count(1)));
        assert_eq!(number(&key::MATE_BREAKEND_IDS), Some(Number::Unknown));
        assert_eq!(number(&key::PARTNER_BREAKEND_ID), Some(Number::Count(1)));
        assert_eq!(number(&key::BREAKEND_EVENT_ID), Some(Number::Count(1)));
        assert_eq!(
            number(&key::BREAKEND_CONFIDENCE_INTERVALS),
            Some(Number::Count(2))
        );
        // assert_eq!(number(&Key::BreakendReadDepth), Some(Number::Count(1)));
        assert_eq!(number(&key::ADJACENT_READ_DEPTHS), Some(Number::Unknown));
        assert_eq!(number(&key::BREAKEND_COPY_NUMBER), Some(Number::Count(1)));
        assert_eq!(number(&key::ADJACENT_COPY_NUMBER), Some(Number::Unknown));
        assert_eq!(
            number(&key::COPY_NUMBER_CONFIDENCE_INTERVALS),
            Some(Number::Count(2))
        );
        assert_eq!(
            number(&key::ADJACENT_COPY_NUMBER_CONFIDENCE_INTERVALS),
            Some(Number::Unknown)
        );

        assert!(number(&"NDLS".parse()?).is_none());

        Ok(())
    }

    #[test]
    fn test_ty() -> Result<(), key::ParseError> {
        assert_eq!(ty(&key::ANCESTRAL_ALLELE), Some(Type::String));
        assert_eq!(ty(&key::ALLELE_COUNT), Some(Type::Integer));
        assert_eq!(ty(&key::TOTAL_READ_DEPTHS), Some(Type::Integer));
        assert_eq!(ty(&key::FORWARD_STRAND_READ_DEPTHS), Some(Type::Integer));
        assert_eq!(ty(&key::REVERSE_STRAND_READ_DEPTHS), Some(Type::Integer));
        assert_eq!(ty(&key::ALLELE_FREQUENCIES), Some(Type::Float));
        assert_eq!(ty(&key::TOTAL_ALLELE_COUNT), Some(Type::Integer));
        assert_eq!(ty(&key::BASE_QUALITY), Some(Type::Float));
        assert_eq!(ty(&key::CIGAR), Some(Type::String));
        assert_eq!(ty(&key::IS_IN_DB_SNP), Some(Type::Flag));
        assert_eq!(ty(&key::TOTAL_DEPTH), Some(Type::Integer));
        // assert_eq!(ty(&END_POSITION), Some(Type::Integer));
        assert_eq!(ty(&key::IS_IN_HAP_MAP_2), Some(Type::Flag));
        assert_eq!(ty(&key::IS_IN_HAP_MAP_3), Some(Type::Flag));
        assert_eq!(ty(&key::MAPPING_QUALITY), Some(Type::Float));
        assert_eq!(ty(&key::ZERO_MAPPING_QUALITY_COUNT), Some(Type::Integer));
        assert_eq!(ty(&key::SAMPLES_WITH_DATA_COUNT), Some(Type::Integer));
        assert_eq!(ty(&key::STRAND_BIAS), Some(Type::Integer));
        assert_eq!(ty(&key::IS_SOMATIC_MUTATION), Some(Type::Flag));
        assert_eq!(ty(&key::IS_VALIDATED), Some(Type::Flag));
        assert_eq!(ty(&key::IS_IN_1000_GENOMES), Some(Type::Flag));

        assert_eq!(ty(&key::IS_IMPRECISE), Some(Type::Flag));
        assert_eq!(ty(&key::IS_NOVEL), Some(Type::Flag));
        assert_eq!(ty(&key::END_POSITION), Some(Type::Integer));
        assert_eq!(ty(&key::SV_TYPE), Some(Type::String));
        assert_eq!(ty(&key::SV_LENGTHS), Some(Type::Integer));
        assert_eq!(ty(&key::POSITION_CONFIDENCE_INTERVALS), Some(Type::Integer));
        assert_eq!(ty(&key::END_CONFIDENCE_INTERVALS), Some(Type::Integer));
        assert_eq!(ty(&key::MICROHOMOLOGY_LENGTHS), Some(Type::Integer));
        assert_eq!(ty(&key::MICROHOMOLOGY_SEQUENCES), Some(Type::String));
        assert_eq!(ty(&key::BREAKPOINT_IDS), Some(Type::String));
        assert_eq!(ty(&key::MOBILE_ELEMENT_INFO), Some(Type::String));
        assert_eq!(
            ty(&key::MOBILE_ELEMENT_TRANSDUCTION_INFO),
            Some(Type::String)
        );
        assert_eq!(ty(&key::DBV_ID), Some(Type::String));
        assert_eq!(ty(&key::DB_VAR_ID), Some(Type::String));
        assert_eq!(ty(&key::DB_RIP_ID), Some(Type::String));
        assert_eq!(ty(&key::MATE_BREAKEND_IDS), Some(Type::String));
        assert_eq!(ty(&key::PARTNER_BREAKEND_ID), Some(Type::String));
        assert_eq!(ty(&key::BREAKEND_EVENT_ID), Some(Type::String));
        assert_eq!(ty(&key::BREAKEND_CONFIDENCE_INTERVALS), Some(Type::Integer));
        // assert_eq!(ty(&Key::BreakendReadDepth), Some(Type::Integer));
        assert_eq!(ty(&key::ADJACENT_READ_DEPTHS), Some(Type::Integer));
        assert_eq!(ty(&key::BREAKEND_COPY_NUMBER), Some(Type::Integer));
        assert_eq!(ty(&key::ADJACENT_COPY_NUMBER), Some(Type::Integer));
        assert_eq!(
            ty(&key::COPY_NUMBER_CONFIDENCE_INTERVALS),
            Some(Type::Integer)
        );
        assert_eq!(
            ty(&key::ADJACENT_COPY_NUMBER_CONFIDENCE_INTERVALS),
            Some(Type::Integer)
        );

        assert!(ty(&"NDLS".parse()?).is_none());

        Ok(())
    }

    #[test]
    fn test_description() -> Result<(), key::ParseError> {
        assert_eq!(
            description(&key::ANCESTRAL_ALLELE),
            Some("Ancestral allele")
        );
        assert_eq!(
            description(&key::ALLELE_COUNT),
            Some("Allele count in genotypes, for each ALT allele, in the same order as listed")
        );
        assert_eq!(
            description(&key::TOTAL_READ_DEPTHS),
            Some("Total read depth for each allele")
        );
        assert_eq!(
            description(&key::FORWARD_STRAND_READ_DEPTHS),
            Some("Read depth for each allele on the forward strand")
        );
        assert_eq!(
            description(&key::REVERSE_STRAND_READ_DEPTHS),
            Some("Read depth for each allele on the reverse strand")
        );
        assert_eq!(
            description(&key::ALLELE_FREQUENCIES),
            Some("Allele frequency for each ALT allele in the same order as listed")
        );
        assert_eq!(
            description(&key::TOTAL_ALLELE_COUNT),
            Some("Total number of alleles in called genotypes")
        );
        assert_eq!(description(&key::BASE_QUALITY), Some("RMS base quality"));
        assert_eq!(
            description(&key::CIGAR),
            Some(
                "Cigar string describing how to align an alternate allele to the reference allele"
            )
        );
        assert_eq!(description(&key::IS_IN_DB_SNP), Some("dbSNP membership"));
        assert_eq!(
            description(&key::TOTAL_DEPTH),
            Some("Combined depth across samples")
        );
        // Self::EndPosition.description(), Some("End position on CHROM"));
        assert_eq!(
            description(&key::IS_IN_HAP_MAP_2),
            Some("HapMap2 membership")
        );
        assert_eq!(
            description(&key::IS_IN_HAP_MAP_3),
            Some("HapMap3 membership")
        );
        assert_eq!(
            description(&key::MAPPING_QUALITY),
            Some("RMS mapping quality")
        );
        assert_eq!(
            description(&key::ZERO_MAPPING_QUALITY_COUNT),
            Some("Number of MAPQ == 0 reads")
        );
        assert_eq!(
            description(&key::SAMPLES_WITH_DATA_COUNT),
            Some("Number of samples with data")
        );
        assert_eq!(description(&key::STRAND_BIAS), Some("Strand bias"));
        assert_eq!(
            description(&key::IS_SOMATIC_MUTATION),
            Some("Somatic mutation")
        );
        assert_eq!(
            description(&key::IS_VALIDATED),
            Some("Validated by follow-up experiment")
        );
        assert_eq!(
            description(&key::IS_IN_1000_GENOMES),
            Some("1000 Genomes membership")
        );

        assert_eq!(
            description(&key::IS_IMPRECISE),
            Some("Imprecise structural variation")
        );
        assert_eq!(
            description(&key::IS_NOVEL),
            Some("Indicates a novel structural variation")
        );
        assert_eq!(
            description(&key::END_POSITION),
            Some("End position of the variant described in this record")
        );
        assert_eq!(
            description(&key::SV_TYPE),
            Some("Type of structural variant")
        );
        assert_eq!(
            description(&key::SV_LENGTHS),
            Some("Difference in length between REF and ALT alleles")
        );
        assert_eq!(
            description(&key::POSITION_CONFIDENCE_INTERVALS),
            Some("Confidence interval around POS for imprecise variants")
        );
        assert_eq!(
            description(&key::END_CONFIDENCE_INTERVALS),
            Some("Confidence interval around END for imprecise variants")
        );
        assert_eq!(
            description(&key::MICROHOMOLOGY_LENGTHS),
            Some("Length of base pair identical micro-homology at event breakpoints")
        );
        assert_eq!(
            description(&key::MICROHOMOLOGY_SEQUENCES),
            Some("Sequence of base pair identical micro-homology at event breakpoints")
        );
        assert_eq!(
            description(&key::BREAKPOINT_IDS),
            Some("ID of the assembled alternate allele in the assembly file")
        );
        assert_eq!(
            description(&key::MOBILE_ELEMENT_INFO),
            Some("Mobile element info of the form NAME,START,END,POLARITY")
        );
        assert_eq!(
            description(&key::MOBILE_ELEMENT_TRANSDUCTION_INFO),
            Some("Mobile element transduction info of the form CHR,START,END,POLARITY")
        );
        assert_eq!(
            description(&key::DBV_ID),
            Some("ID of this element in Database of Genomic Variation")
        );
        assert_eq!(
            description(&key::DB_VAR_ID),
            Some("ID of this element in DBVAR")
        );
        assert_eq!(
            description(&key::DB_RIP_ID),
            Some("ID of this element in DBRIP")
        );
        assert_eq!(
            description(&key::MATE_BREAKEND_IDS),
            Some("ID of mate breakends")
        );
        assert_eq!(
            description(&key::PARTNER_BREAKEND_ID),
            Some("ID of partner breakend")
        );
        assert_eq!(
            description(&key::BREAKEND_EVENT_ID),
            Some("ID of event associated to breakend")
        );
        assert_eq!(
            description(&key::BREAKEND_CONFIDENCE_INTERVALS),
            Some("Confidence interval around the inserted material between breakends")
        );
        // Self::BreakendReadDepth.description(), Some("Read Depth of segment containing breakend"));
        assert_eq!(
            description(&key::ADJACENT_READ_DEPTHS),
            Some("Read Depth of adjacency")
        );
        assert_eq!(
            description(&key::BREAKEND_COPY_NUMBER),
            Some("Copy number of segment containing breakend")
        );
        assert_eq!(
            description(&key::ADJACENT_COPY_NUMBER),
            Some("Copy number of adjacency")
        );
        assert_eq!(
            description(&key::COPY_NUMBER_CONFIDENCE_INTERVALS),
            Some("Confidence interval around copy number for the segment")
        );
        assert_eq!(
            description(&key::ADJACENT_COPY_NUMBER_CONFIDENCE_INTERVALS),
            Some("Confidence interval around copy number for the adjacency")
        );

        assert!(description(&"NDLS".parse()?).is_none());

        Ok(())
    }
}
