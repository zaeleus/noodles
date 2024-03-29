//! Lazily-evaluated BCF record and fields.

mod convert;
mod filters;
mod genotypes;
mod info;
pub(crate) mod value;

pub(crate) use self::value::Value;
pub use self::{filters::Filters, genotypes::Genotypes, info::Info};

use std::io;

use noodles_vcf as vcf;

/// A chromosome ID.
pub type ChromosomeId = usize;

/// A BCF record.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    pub(crate) chrom: ChromosomeId,
    pub(crate) pos: vcf::record::Position,
    pub(crate) rlen: usize,
    pub(crate) qual: Option<vcf::record::QualityScore>,
    pub(crate) id: vcf::record::Ids,
    pub(crate) r#ref: vcf::record::ReferenceBases,
    pub(crate) alt: vcf::record::AlternateBases,
    pub(crate) filter: Filters,
    pub(crate) info: Info,
    pub(crate) genotypes: Genotypes,
}

impl Record {
    /// Returns the chromosome ID of the record.
    ///
    /// The chromosome ID represents an index in the contig string map, which associates an ID (by
    /// position) with a contig record in the VCF header (by name). That is, to get the associated
    /// contig record in the VCF header, the contig string map must first be queried by position to
    /// find the chromosome name, and then the contigs in the VCF header can be queried by name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let record = bcf::lazy::Record::default();
    /// assert_eq!(record.chromosome_id(), 0);
    /// ```
    pub fn chromosome_id(&self) -> ChromosomeId {
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
    /// let record = bcf::lazy::Record::default();
    /// assert_eq!(usize::from(record.position()), 1);
    /// ```
    pub fn position(&self) -> vcf::record::Position {
        self.pos
    }

    pub(crate) fn rlen(&self) -> usize {
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
    /// let record = bcf::lazy::Record::default();
    /// assert_eq!(record.end().map(usize::from)?, 1);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn end(&self) -> io::Result<vcf::record::Position> {
        use vcf::record::Position;

        let start = usize::from(self.position());
        let len = self.rlen();
        let end = start + len - 1;

        Ok(Position::from(end))
    }

    /// Returns the quality score.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let record = bcf::lazy::Record::default();
    /// assert!(record.quality_score().is_none());
    /// ```
    pub fn quality_score(&self) -> Option<vcf::record::QualityScore> {
        self.qual
    }

    /// Returns the IDs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let record = bcf::lazy::Record::default();
    /// assert!(record.ids().is_empty());
    /// ```
    pub fn ids(&self) -> &vcf::record::Ids {
        &self.id
    }

    pub(crate) fn reference_bases(&self) -> &vcf::record::ReferenceBases {
        &self.r#ref
    }

    pub(crate) fn alternate_bases(&self) -> &vcf::record::AlternateBases {
        &self.alt
    }

    /// Returns the filters.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let record = bcf::lazy::Record::default();
    /// assert!(record.filters().is_empty());
    /// ```
    pub fn filters(&self) -> &Filters {
        &self.filter
    }

    /// Returns the info.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let record = bcf::lazy::Record::default();
    /// assert!(record.info().is_empty());
    /// ```
    pub fn info(&self) -> &Info {
        &self.info
    }

    /// Returns the genotypes.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let record = bcf::lazy::Record::default();
    /// assert!(record.genotypes().is_empty());
    /// ```
    pub fn genotypes(&self) -> &Genotypes {
        &self.genotypes
    }
}

impl Default for Record {
    fn default() -> Self {
        use vcf::record::reference_bases::Base;

        Self {
            chrom: 0,
            pos: vcf::record::Position::from(1),
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
