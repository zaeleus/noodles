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
    chrom: ChromosomeId,
    pos: vcf::record::Position,
    rlen: usize,
    qual: Option<vcf::record::QualityScore>,
    id: vcf::record::Ids,
    pub(crate) r#ref: vcf::record::ReferenceBases,
    pub(crate) alt: vcf::record::AlternateBases,
    filter: Filters,
    info: Info,
    genotypes: Genotypes,
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

    pub(crate) fn chromosome_id_mut(&mut self) -> &mut ChromosomeId {
        &mut self.chrom
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

    /// Returns a mutable reference to the start position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// use noodles_vcf::record::Position;
    ///
    /// let mut record = bcf::lazy::Record::default();
    /// *record.position_mut() = Position::from(8);
    ///
    /// assert_eq!(usize::from(record.position()), 8);
    /// ```
    pub fn position_mut(&mut self) -> &mut vcf::record::Position {
        &mut self.pos
    }

    pub(crate) fn rlen(&self) -> usize {
        self.rlen
    }

    pub(crate) fn rlen_mut(&mut self) -> &mut usize {
        &mut self.rlen
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

    /// Return a mutable reference to the quality score.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// use noodles_vcf::record::QualityScore;
    ///
    /// let mut record = bcf::lazy::Record::default();
    /// *record.quality_score_mut() = QualityScore::try_from(13.0).map(Some)?;
    ///
    /// assert_eq!(record.quality_score().map(f32::from), Some(13.0));
    /// # Ok::<_, noodles_vcf::record::quality_score::TryFromFloatError>(())
    /// ```
    pub fn quality_score_mut(&mut self) -> &mut Option<vcf::record::QualityScore> {
        &mut self.qual
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

    /// Returns a mutable reference to the IDs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// use noodles_vcf::record::Ids;
    ///
    /// let mut record = bcf::lazy::Record::default();
    /// let ids: Ids = "nd0".parse()?;
    /// *record.ids_mut() = ids.clone();
    ///
    /// assert_eq!(record.ids(), &ids);
    /// # Ok::<_, noodles_vcf::record::ids::ParseError>(())
    /// ```
    pub fn ids_mut(&mut self) -> &mut vcf::record::Ids {
        &mut self.id
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

    pub(crate) fn filters_mut(&mut self) -> &mut Filters {
        &mut self.filter
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

    pub(crate) fn info_mut(&mut self) -> &mut Info {
        &mut self.info
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

    pub(crate) fn genotypes_mut(&mut self) -> &mut Genotypes {
        &mut self.genotypes
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
