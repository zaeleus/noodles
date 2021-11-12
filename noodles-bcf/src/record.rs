//! BCF record and fields.

mod convert;
mod filters;
mod genotypes;
mod info;
pub mod value;

pub use self::{filters::Filters, genotypes::Genotypes, info::Info, value::Value};

use std::io;

use noodles_vcf as vcf;

/// A BCF record.
///
/// A `bcf::Record` wraps a raw byte buffer, and the fields should be considered immutable.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    pub(crate) chrom: i32,
    pub(crate) pos: vcf::record::Position,
    pub(crate) rlen: i32,
    pub(crate) qual: Option<vcf::record::QualityScore>,
    pub(crate) id: vcf::record::Ids,
    pub(crate) r#ref: vcf::record::ReferenceBases,
    pub(crate) alt: vcf::record::AlternateBases,
    filter: Filters,
    info: Info,
    genotypes: Genotypes,
}

impl Record {
    /// Returns the chromosome ID of the record.
    ///
    /// The chromosome ID is the index of the associated contig in the VCF header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let record = bcf::Record::default();
    /// assert_eq!(record.chromosome_id(), 0);
    /// ```
    pub fn chromosome_id(&self) -> i32 {
        self.chrom
    }

    /// Returns the start position of this record.
    ///
    /// Despite the BCF format using 0-based positions, this normalizes the value as a 1-based
    /// position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let record = bcf::Record::default();
    /// assert_eq!(i32::from(record.position()), 1);
    /// ```
    pub fn position(&self) -> vcf::record::Position {
        self.pos
    }

    fn rlen(&self) -> i32 {
        self.rlen
    }

    /// Returns the end position of this record.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// let record = bcf::Record::default();
    /// assert_eq!(record.end().map(i32::from)?, 1);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn end(&self) -> io::Result<vcf::record::Position> {
        use vcf::record::Position;

        let start = i32::from(self.position());
        let len = self.rlen();
        let end = start + len - 1;

        Position::try_from(end).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    pub(crate) fn quality_score(&self) -> Option<vcf::record::QualityScore> {
        self.qual
    }

    pub(crate) fn ids(&self) -> &vcf::record::Ids {
        &self.id
    }

    pub(crate) fn reference_bases(&self) -> &vcf::record::ReferenceBases {
        &self.r#ref
    }

    pub(crate) fn alternate_bases(&self) -> &vcf::record::AlternateBases {
        &self.alt
    }

    pub(crate) fn filters(&self) -> &Filters {
        &self.filter
    }

    pub(crate) fn filters_mut(&mut self) -> &mut Filters {
        &mut self.filter
    }

    pub(crate) fn info(&self) -> &Info {
        &self.info
    }

    pub(crate) fn info_mut(&mut self) -> &mut Info {
        &mut self.info
    }

    pub(crate) fn genotypes(&self) -> &Genotypes {
        &self.genotypes
    }

    pub(crate) fn genotypes_mut(&mut self) -> &mut Genotypes {
        &mut self.genotypes
    }
}

impl Default for Record {
    fn default() -> Self {
        use vcf::record::reference_bases::Base;

        Self {
            chrom: 0,
            pos: vcf::record::Position::try_from(1).unwrap(),
            rlen: 1,
            qual: None,
            id: vcf::record::Ids::default(),
            r#ref: vcf::record::ReferenceBases::try_from(vec![Base::A]).unwrap(),
            alt: vcf::record::AlternateBases::default(),
            filter: Filters::default(),
            info: Info::default(),
            genotypes: Genotypes::default(),
        }
    }
}
