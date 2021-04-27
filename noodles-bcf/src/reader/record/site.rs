use noodles_vcf::record::{Filters, Ids, Info, QualityScore};

#[derive(Clone, Debug, PartialEq)]
pub struct Site {
    pub chrom: i32,
    pub pos: i32,
    pub rlen: i32,
    pub qual: QualityScore,
    pub n_allele_info: i32,
    pub n_fmt_sample: u32,
    pub id: Ids,
    pub ref_alt: Vec<String>,
    pub filter: Filters,
    pub info: Info,
}
