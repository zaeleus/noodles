//! BCF record.

pub(crate) mod codec;
mod convert;
mod fields;
mod filters;
mod genotypes;
mod ids;
mod info;
mod reference_bases;
mod value;

use std::io;

use noodles_vcf as vcf;

use self::fields::Fields;
pub(crate) use self::value::Value;
pub use self::{
    filters::Filters, genotypes::Genotypes, ids::Ids, info::Info, reference_bases::ReferenceBases,
};

/// A chromosome ID.
pub type ChromosomeId = usize;

/// A BCF record.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    fields: Fields,
    pub(crate) chrom: ChromosomeId,
    pub(crate) pos: vcf::record::Position,
    pub(crate) rlen: usize,
    pub(crate) qual: Option<f32>,
    pub(crate) id: vcf::record::Ids,
    pub(crate) r#ref: String,
    pub(crate) alt: vcf::record::AlternateBases,
    pub(crate) filter: Filters,
    pub(crate) info: Info,
}

impl Record {
    pub(crate) fn fields_mut(&mut self) -> &mut Fields {
        &mut self.fields
    }

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
    /// let record = bcf::Record::default();
    /// assert_eq!(record.chromosome_id()?, 0);
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn chromosome_id(&self) -> io::Result<usize> {
        let n = self.fields.reference_sequence_id();
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
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
    /// assert_eq!(record.position().map(usize::from)?, 1);
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn position(&self) -> io::Result<vcf::record::Position> {
        let n = self.fields.position();

        usize::try_from(n)
            .map(|m| m + 1)
            .map(vcf::record::Position::from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    pub(crate) fn rlen(&self) -> io::Result<usize> {
        let n = self.fields.span();
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
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
    /// assert_eq!(record.end().map(usize::from)?, 1);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn end(&self) -> io::Result<vcf::record::Position> {
        use vcf::record::Position;

        let start = self.position().map(usize::from)?;
        let len = self.rlen()?;
        let end = start + len - 1;

        Ok(Position::from(end))
    }

    /// Returns the quality score.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let record = bcf::Record::default();
    /// assert!(record.quality_score()?.is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn quality_score(&self) -> io::Result<Option<f32>> {
        self.fields.quality_score()
    }

    /// Returns the IDs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let record = bcf::Record::default();
    /// assert!(record.ids().is_empty());
    /// ```
    pub fn ids(&self) -> Ids<'_> {
        self.fields.ids()
    }

    pub(crate) fn reference_bases(&self) -> ReferenceBases<'_> {
        self.fields.reference_bases()
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
    /// let record = bcf::Record::default();
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
    /// let record = bcf::Record::default();
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
    /// let record = bcf::Record::default();
    /// assert!(record.genotypes()?.is_empty());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn genotypes(&self) -> io::Result<Genotypes<'_>> {
        self.fields.genotypes()
    }
}

impl Default for Record {
    fn default() -> Self {
        Self {
            fields: Fields::default(),
            chrom: 0,
            pos: vcf::record::Position::from(1),
            rlen: 1,
            qual: None,
            id: vcf::record::Ids::default(),
            r#ref: String::from("N"),
            alt: vcf::record::AlternateBases::default(),
            filter: Filters::default(),
            info: Info::default(),
        }
    }
}
