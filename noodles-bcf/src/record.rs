//! BCF record.

mod alternate_bases;
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
    alternate_bases::AlternateBases, filters::Filters, genotypes::Genotypes, ids::Ids, info::Info,
    reference_bases::ReferenceBases,
};

/// A chromosome ID.
pub type ChromosomeId = usize;

/// A BCF record.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Record(Fields);

impl Record {
    pub(crate) fn fields_mut(&mut self) -> &mut Fields {
        &mut self.0
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
        let n = self.0.reference_sequence_id();
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
        let n = self.0.position();

        usize::try_from(n)
            .map(|m| m + 1)
            .map(vcf::record::Position::from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    pub(crate) fn rlen(&self) -> io::Result<usize> {
        let n = self.0.span();
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
        self.0.quality_score()
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
        self.0.ids()
    }

    /// Returns the reference bases.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let record = bcf::Record::default();
    /// assert_eq!(record.reference_bases().as_ref(), b"N");
    /// ```
    pub fn reference_bases(&self) -> ReferenceBases<'_> {
        self.0.reference_bases()
    }

    /// Returns the alternate bases.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let record = bcf::Record::default();
    /// assert!(record.alternate_bases().is_empty());
    /// ```
    pub fn alternate_bases(&self) -> AlternateBases<'_> {
        self.0.alternate_bases()
    }

    /// Returns the filters.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let record = bcf::Record::default();
    /// assert!(record.filters().is_empty()?);
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn filters(&self) -> Filters<'_> {
        self.0.filters()
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
    pub fn info(&self) -> Info<'_> {
        self.0.info()
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
        self.0.genotypes()
    }
}
