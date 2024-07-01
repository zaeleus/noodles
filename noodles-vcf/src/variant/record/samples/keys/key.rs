//! Variant record samples key.

/// Read depth for each allele (`AD`).
pub const READ_DEPTHS: &str = "AD";

/// Read depth for each allele on the forward strand (`ADF`).
pub const FORWARD_STRAND_READ_DEPTHS: &str = "ADF";

/// Read depth for each allele on the reverse strand (`ADR`).
pub const REVERSE_STRAND_READ_DEPTHS: &str = "ADR";

/// Read depth (`DP`).
pub const READ_DEPTH: &str = "DP";

/// Expected alternate allele counts (`EC`).
pub const EXPECTED_ALTERNATE_ALLELE_COUNTS: &str = "EC";

/// Length of <*> reference block (`LEN`).
///
/// Added in VCF 4.5.
pub const LENGTH: &str = "LEN";

/// Filter indicating if this genotype was "called" (`FT`).
pub const FILTER: &str = "FT";

/// Genotype likelihoods (`GL`).
pub const GENOTYPE_LIKELIHOODS: &str = "GL";

/// Genotype posterior probabilities (`GP`).
pub const GENOTYPE_POSTERIOR_PROBABILITIES: &str = "GP";

/// Conditional genotype quality (`GQ`).
pub const CONDITIONAL_GENOTYPE_QUALITY: &str = "GQ";

/// Genotype (`GT`).
pub const GENOTYPE: &str = "GT";

/// Haplotype quality (`HQ`).
pub const HAPLOTYPE_QUALITY: &str = "HQ";

/// Reserved (`LA`).
///
/// Added in VCF 4.5.
pub const RESERVED_LA: &str = "LA";

/// 1-based indices into ALT, indicating which alleles are relevant (local) for the current sample
/// (`LAA`).
///
/// Added in VCF 4.5.
pub const LOCAL_ALTERNATIVE_ALLELE: &str = "LAA";

/// Local-allele representation of AD (`LAD`).
///
/// Added in VCF 4.5.
pub const LOCAL_READ_DEPTHS: &str = "LAD";

/// Local-allele representation of ADF (`LADF`).
///
/// Added in VCF 4.5.
pub const LOCAL_FORWARD_STRAND_READ_DEPTHS: &str = "LADF";

/// Local-allele representation of ADR (`LADR`).
///
/// Added in VCF 4.5.
pub const LOCAL_REVERSE_STRAND_READ_DEPTHS: &str = "LADR";

/// Local-allele representation of EC (`LEC`).
///
/// Added in VCF 4.5.
pub const LOCAL_EXPECTED_ALTERNATE_ALLELE_COUNTS: &str = "LEC";

/// Local-allele representation of GL (`LGL`).
///
/// Added in VCF 4.5.
pub const LOCAL_GENOTYPE_LIKELIHOODS: &str = "LGL";

/// Local-allele representation of GP (`LGP`).
///
/// Added in VCF 4.5.
pub const LOCAL_GENOTYPE_POSTERIOR_PROBABILITIES: &str = "LGP";

/// Local-allele representation of PL (`LPL`).
///
/// Added in VCF 4.5.
pub const LOCAL_ROUNDED_GENOTYPE_LIKELIHOODS: &str = "LPL";

/// Local-allele representation of PP (`LPP`).
///
/// Added in VCF 4.5.
pub const LOCAL_ROUNDED_GENOTYPE_POSTERIOR_PROBABILITIES: &str = "LPP";

/// RMS mapping quality (`MQ`).
pub const MAPPING_QUALITY: &str = "MQ";

/// Phred-scaled genotype likelihoods rounded to the closest integer (`PL`).
pub const ROUNDED_GENOTYPE_LIKELIHOODS: &str = "PL";

/// Phred-scaled genotype posterior probabilities rounded to the closest integer (`PP`).
pub const ROUNDED_GENOTYPE_POSTERIOR_PROBABILITIES: &str = "PP";

/// Phasing quality (`PQ`).
pub const PHASING_QUALITY: &str = "PQ";

/// Phase set (`PS`).
pub const PHASE_SET: &str = "PS";

/// Phase set list (`PSL`).
///
/// Added in VCF 4.4.
pub const PHASE_SET_LIST: &str = "PSL";

/// Phase set list ordinal (`PSO`).
///
/// Added in VCF 4.4.
pub const PHASE_SET_LIST_ORDINALS: &str = "PSO";

/// Phase set list quality (`PSQ`).
///
/// Added in VCF 4.4.
pub const PHASE_SET_LIST_QUALITIES: &str = "PSQ";

/// Copy number genotype for imprecise events (`CN`).
pub const GENOTYPE_COPY_NUMBER: &str = "CN";

/// Confidence interval around copy number (`CICN`).
///
/// Added in VCF 4.4.
pub const COPY_NUMBER_CONFIDENCE_INTERVAL: &str = "CICN";

/// Copy number genotype quality for imprecise events (`CNQ`).
pub const GENOTYPE_COPY_NUMBER_QUALITY: &str = "CNQ";

/// Copy number genotype likelihood for imprecise events (`CNL`).
pub const GENOTYPE_COPY_NUMBER_LIKELIHOODS: &str = "CNL";

/// Copy number posterior probabilities (`CNP`).
pub const GENOTYPE_COPY_NUMBER_POSTERIOR_PROBABILITIES: &str = "CNP";

/// Phred style probability score that the variant is novel (`NQ`).
pub const NOVEL_VARIANT_QUALITY_SCORE: &str = "NQ";

/// Unique haplotype identifier (`HAP`).
pub const HAPLOTYPE_ID: &str = "HAP";

/// Unique identifier of ancestral haplotype (`AHAP`).
pub const ANCESTRAL_HAPLOTYPE_ID: &str = "AHAP";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(READ_DEPTHS, "AD");
        assert_eq!(FORWARD_STRAND_READ_DEPTHS, "ADF");
        assert_eq!(REVERSE_STRAND_READ_DEPTHS, "ADR");
        assert_eq!(READ_DEPTH, "DP");
        assert_eq!(EXPECTED_ALTERNATE_ALLELE_COUNTS, "EC");
        assert_eq!(LENGTH, "LEN");
        assert_eq!(FILTER, "FT");
        assert_eq!(GENOTYPE_LIKELIHOODS, "GL");
        assert_eq!(GENOTYPE_POSTERIOR_PROBABILITIES, "GP");
        assert_eq!(CONDITIONAL_GENOTYPE_QUALITY, "GQ");
        assert_eq!(GENOTYPE, "GT");
        assert_eq!(HAPLOTYPE_QUALITY, "HQ");
        assert_eq!(RESERVED_LA, "LA");
        assert_eq!(LOCAL_ALTERNATIVE_ALLELE, "LAA");
        assert_eq!(LOCAL_READ_DEPTHS, "LAD");
        assert_eq!(LOCAL_FORWARD_STRAND_READ_DEPTHS, "LADF");
        assert_eq!(LOCAL_REVERSE_STRAND_READ_DEPTHS, "LADR");
        assert_eq!(LOCAL_EXPECTED_ALTERNATE_ALLELE_COUNTS, "LEC");
        assert_eq!(LOCAL_GENOTYPE_LIKELIHOODS, "LGL");
        assert_eq!(LOCAL_GENOTYPE_POSTERIOR_PROBABILITIES, "LGP");
        assert_eq!(LOCAL_ROUNDED_GENOTYPE_LIKELIHOODS, "LPL");
        assert_eq!(LOCAL_ROUNDED_GENOTYPE_POSTERIOR_PROBABILITIES, "LPP");
        assert_eq!(MAPPING_QUALITY, "MQ");
        assert_eq!(ROUNDED_GENOTYPE_LIKELIHOODS, "PL");
        assert_eq!(ROUNDED_GENOTYPE_POSTERIOR_PROBABILITIES, "PP");
        assert_eq!(PHASING_QUALITY, "PQ");
        assert_eq!(PHASE_SET, "PS");
        assert_eq!(PHASE_SET_LIST, "PSL");
        assert_eq!(PHASE_SET_LIST_ORDINALS, "PSO");
        assert_eq!(PHASE_SET_LIST_QUALITIES, "PSQ");

        assert_eq!(GENOTYPE_COPY_NUMBER, "CN");
        assert_eq!(COPY_NUMBER_CONFIDENCE_INTERVAL, "CICN");
        assert_eq!(GENOTYPE_COPY_NUMBER_QUALITY, "CNQ");
        assert_eq!(GENOTYPE_COPY_NUMBER_LIKELIHOODS, "CNL");
        assert_eq!(GENOTYPE_COPY_NUMBER_POSTERIOR_PROBABILITIES, "CNP");
        assert_eq!(NOVEL_VARIANT_QUALITY_SCORE, "NQ");
        assert_eq!(HAPLOTYPE_ID, "HAP");
        assert_eq!(ANCESTRAL_HAPLOTYPE_ID, "AHAP");
    }
}
