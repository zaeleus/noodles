//! VCF record info key.

/// Ancestral allele (`AA`).
pub const ANCESTRAL_ALLELE: &str = "AA";

/// Allele count in genotypes, for each ALT allele, in the same order as listed (`AC`).
pub const ALLELE_COUNT: &str = "AC";

/// Total read depth for each allele (`AD`).
pub const TOTAL_READ_DEPTHS: &str = "AD";

/// Read depth for each allele on the forward strand (`ADF`).
pub const FORWARD_STRAND_READ_DEPTHS: &str = "ADF";

/// Read depth for each allele on the reverse strand (`ADR`).
pub const REVERSE_STRAND_READ_DEPTHS: &str = "ADR";

/// Allele frequency for each ALT allele in the same order as listed (`AF`).
pub const ALLELE_FREQUENCIES: &str = "AF";

/// Total number of alleles in called genotypes (`AN`).
pub const TOTAL_ALLELE_COUNT: &str = "AN";

/// RMS base quality (`BQ`).
pub const BASE_QUALITY: &str = "BQ";

/// Cigar string describing how to align an alternate allele to the reference allele (`CIGAR`).
pub const CIGAR: &str = "CIGAR";

/// dbSNP membership (`DB`).
pub const IS_IN_DB_SNP: &str = "DB";

/// Combined depth across samples (`DP`).
pub const TOTAL_DEPTH: &str = "DP";

/// HapMap2 membership (`H2`).
pub const IS_IN_HAP_MAP_2: &str = "H2";

/// HapMap3 membership (`H3`).
pub const IS_IN_HAP_MAP_3: &str = "H3";

/// RMS mapping quality (`MQ`).
pub const MAPPING_QUALITY: &str = "MQ";

/// Number of MAPQ == 0 reads (`MQ0`).
pub const ZERO_MAPPING_QUALITY_COUNT: &str = "MQ0";

/// Number of samples with data (`NS`).
pub const SAMPLES_WITH_DATA_COUNT: &str = "NS";

/// Strand bias (`SB`).
pub const STRAND_BIAS: &str = "SB";

/// Somatic mutation (`SOMATIC`).
pub const IS_SOMATIC_MUTATION: &str = "SOMATIC";

/// Validated by follow-up experiment (`VALIDATED`).
pub const IS_VALIDATED: &str = "VALIDATED";

/// 1000 Genomes membership (`1000G`).
pub const IS_IN_1000_GENOMES: &str = "1000G";

/// Imprecise structural variation (`IMPRECISE`).
pub const IS_IMPRECISE: &str = "IMPRECISE";

/// Indicates a novel structural variation (`NOVEL`).
pub const IS_NOVEL: &str = "NOVEL";

/// End position of the variant described in this record (`END`).
pub const END_POSITION: &str = "END";

/// Type of structural variant (`SVTYPE`).
///
/// Deprecated in VCF 4.4.
pub const SV_TYPE: &str = "SVTYPE";

/// Difference in length between REF and ALT alleles (`SVLEN`).
pub const SV_LENGTHS: &str = "SVLEN";

/// Confidence interval around POS for imprecise variants (`CIPOS`).
pub const POSITION_CONFIDENCE_INTERVALS: &str = "CIPOS";

/// Confidence interval around END for imprecise variants (`CIEND`).
pub const END_CONFIDENCE_INTERVALS: &str = "CIEND";

/// Length of base pair identical micro-homology at event breakpoints (`HOMLEN`).
pub const MICROHOMOLOGY_LENGTHS: &str = "HOMLEN";

/// Sequence of base pair identical micro-homology at event breakpoints (`HOMSEQ`).
pub const MICROHOMOLOGY_SEQUENCES: &str = "HOMSEQ";

/// ID of the assembled alternate allele in the assembly file (`BKPTID`).
pub const BREAKPOINT_IDS: &str = "BKPTID";

/// Mobile element info of the form NAME,START,END,POLARITY (`MEINFO`).
pub const MOBILE_ELEMENT_INFO: &str = "MEINFO";

/// Mobile element transduction info of the form CHR,START,END,POLARITY (`METRANS`).
pub const MOBILE_ELEMENT_TRANSDUCTION_INFO: &str = "METRANS";

/// ID of this element in Database of Genomic Variation (`DBVID`).
pub const DBV_ID: &str = "DBVID";

/// ID of this element in DBVAR (`DBVARID`).
pub const DB_VAR_ID: &str = "DBVARID";

/// ID of this element in DBRIP (`DBRIPID`).
pub const DB_RIP_ID: &str = "DBRIPID";

/// ID of mate breakends (`MATEID`).
pub const MATE_BREAKEND_IDS: &str = "MATEID";

/// ID of partner breakend (`PARID`).
pub const PARTNER_BREAKEND_ID: &str = "PARID";

/// ID of event associated to breakend (`EVENT`).
pub const BREAKEND_EVENT_ID: &str = "EVENT";

/// Type of associated event (`EVENTTYPE`).
///
/// Added in VCF 4.4.
pub const EVENT_TYPE: &str = "EVENTTYPE";

/// Confidence interval around the inserted material between breakends (`CILEN`).
pub const BREAKEND_CONFIDENCE_INTERVALS: &str = "CILEN";

/// Read Depth of adjacency (`DPADJ`).
///
/// Removed in VCF 4.4.
pub const ADJACENT_READ_DEPTHS: &str = "DPADJ";

/// Copy number of segment containing breakend (`CN`).
pub const BREAKEND_COPY_NUMBER: &str = "CN";

/// Copy number of adjacency (`CNADJ`).
///
/// Removed in VCF 4.4.
pub const ADJACENT_COPY_NUMBER: &str = "CNADJ";

/// Confidence interval around copy number for the segment (`CICN`).
pub const COPY_NUMBER_CONFIDENCE_INTERVALS: &str = "CICN";

/// Confidence interval around copy number for the adjacency (`CICNADJ`).
///
/// Removed in VCF 4.4.
pub const ADJACENT_COPY_NUMBER_CONFIDENCE_INTERVALS: &str = "CICNADJ";

/// Claim made by the structural variant call (`SVCLAIM`).
///
/// Valid values are D, J, DJ for abundance, adjacency and both respectively.
///
/// Added in VCF 4.4.
pub const SV_CLAIM: &str = "SVCLAIM";

/// Total number of repeat sequences in this allele (`RN`).
///
/// Added in VCF 4.4.
pub const TOTAL_REPEAT_SEQUENCE_COUNTS: &str = "RN";

/// Repeat unit sequence of the corresponding repeat sequence (`RUS`).
///
/// Added in VCF 4.4.
pub const REPEAT_UNIT_SEQUENCES: &str = "RUS";

/// Repeat unit length of the corresponding repeat sequence (`RUL`).
///
/// Added in VCF 4.4.
pub const REPEAT_UNIT_LENGTHS: &str = "RUL";

/// Repeat unit count of corresponding repeat sequence (`RUC`).
///
/// Added in VCF 4.4.
pub const REPEAT_UNIT_COUNTS: &str = "RUC";

/// Total number of bases in the corresponding repeat sequence (`RB`).
///
/// Added in VCF 4.4.
pub const TOTAL_REPEAT_SEQUENCE_BASE_COUNTS: &str = "RB";

/// Confidence interval around RUC (`CIRUC`).
///
/// Added in VCF 4.4.
pub const REPEAT_UNIT_COUNT_CONFIDENCE_INTERVALS: &str = "CIRUC";

/// Confidence interval around RB (`CIRB`).
///
/// Added in VCF 4.4.
pub const TOTAL_REPEAT_SEQUENCE_BASE_COUNT_CONFIDENCE_INTERVALS: &str = "CIRB";

/// Number of bases in each individual repeat unit (`RUB`).
///
/// Added in VCF 4.4.
pub const REPEAT_UNIT_BASE_COUNTS: &str = "RUB";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_values() {
        assert_eq!(ANCESTRAL_ALLELE, "AA");
        assert_eq!(ALLELE_COUNT, "AC");
        assert_eq!(TOTAL_READ_DEPTHS, "AD");
        assert_eq!(FORWARD_STRAND_READ_DEPTHS, "ADF");
        assert_eq!(REVERSE_STRAND_READ_DEPTHS, "ADR");
        assert_eq!(ALLELE_FREQUENCIES, "AF");
        assert_eq!(TOTAL_ALLELE_COUNT, "AN");
        assert_eq!(BASE_QUALITY, "BQ");
        assert_eq!(CIGAR, "CIGAR");
        assert_eq!(IS_IN_DB_SNP, "DB");
        assert_eq!(TOTAL_DEPTH, "DP");
        assert_eq!(IS_IN_HAP_MAP_2, "H2");
        assert_eq!(IS_IN_HAP_MAP_3, "H3");
        assert_eq!(MAPPING_QUALITY, "MQ");
        assert_eq!(ZERO_MAPPING_QUALITY_COUNT, "MQ0");
        assert_eq!(SAMPLES_WITH_DATA_COUNT, "NS");
        assert_eq!(STRAND_BIAS, "SB");
        assert_eq!(IS_SOMATIC_MUTATION, "SOMATIC");
        assert_eq!(IS_VALIDATED, "VALIDATED");
        assert_eq!(IS_IN_1000_GENOMES, "1000G");

        assert_eq!(IS_IMPRECISE, "IMPRECISE");
        assert_eq!(IS_NOVEL, "NOVEL");
        assert_eq!(END_POSITION, "END");
        assert_eq!(SV_TYPE, "SVTYPE");
        assert_eq!(SV_LENGTHS, "SVLEN");
        assert_eq!(POSITION_CONFIDENCE_INTERVALS, "CIPOS");
        assert_eq!(END_CONFIDENCE_INTERVALS, "CIEND");
        assert_eq!(MICROHOMOLOGY_LENGTHS, "HOMLEN");
        assert_eq!(MICROHOMOLOGY_SEQUENCES, "HOMSEQ");
        assert_eq!(BREAKPOINT_IDS, "BKPTID");
        assert_eq!(MOBILE_ELEMENT_INFO, "MEINFO");
        assert_eq!(MOBILE_ELEMENT_TRANSDUCTION_INFO, "METRANS");
        assert_eq!(DBV_ID, "DBVID");
        assert_eq!(DB_VAR_ID, "DBVARID");
        assert_eq!(DB_RIP_ID, "DBRIPID");
        assert_eq!(MATE_BREAKEND_IDS, "MATEID");
        assert_eq!(PARTNER_BREAKEND_ID, "PARID");
        assert_eq!(BREAKEND_EVENT_ID, "EVENT");
        assert_eq!(BREAKEND_CONFIDENCE_INTERVALS, "CILEN");
        assert_eq!(ADJACENT_READ_DEPTHS, "DPADJ");
        assert_eq!(BREAKEND_COPY_NUMBER, "CN");
        assert_eq!(ADJACENT_COPY_NUMBER, "CNADJ");
        assert_eq!(COPY_NUMBER_CONFIDENCE_INTERVALS, "CICN");
        assert_eq!(ADJACENT_COPY_NUMBER_CONFIDENCE_INTERVALS, "CICNADJ");
        assert_eq!(SV_CLAIM, "SVCLAIM");
        assert_eq!(TOTAL_REPEAT_SEQUENCE_COUNTS, "RN");
        assert_eq!(REPEAT_UNIT_SEQUENCES, "RUS");
        assert_eq!(REPEAT_UNIT_LENGTHS, "RUL");
        assert_eq!(REPEAT_UNIT_COUNTS, "RUC");
        assert_eq!(TOTAL_REPEAT_SEQUENCE_BASE_COUNTS, "RB");
        assert_eq!(REPEAT_UNIT_COUNT_CONFIDENCE_INTERVALS, "CIRUC");
        assert_eq!(
            TOTAL_REPEAT_SEQUENCE_BASE_COUNT_CONFIDENCE_INTERVALS,
            "CIRB"
        );
        assert_eq!(REPEAT_UNIT_BASE_COUNTS, "RUB");
    }
}
