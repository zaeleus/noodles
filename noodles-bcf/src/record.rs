//! BCF record.

mod alternate_bases;
pub(crate) mod codec;
mod fields;
mod filters;
mod ids;
mod info;
mod reference_bases;
pub mod samples;
mod value;

use std::{fmt, io, str};

use noodles_core::Position;
use noodles_vcf::{self as vcf, header::StringMaps};

use self::fields::Fields;
pub(crate) use self::value::Value;
pub use self::{
    alternate_bases::AlternateBases, filters::Filters, ids::Ids, info::Info,
    reference_bases::ReferenceBases, samples::Samples,
};

/// A BCF record.
#[derive(Clone, Default, PartialEq)]
pub struct Record(Fields);

impl Record {
    pub(crate) fn fields_mut(&mut self) -> &mut Fields {
        &mut self.0
    }

    /// Returns the reference sequence ID of the record.
    ///
    /// The reference sequence ID represents an index in the contig string map, which associates an
    /// ID (by position) with a contig record in the VCF header (by name). That is, to get the
    /// associated contig record in the VCF header, the contig string map must first be queried by
    /// position to find the chromosome name, and then the contigs in the VCF header can be queried
    /// by name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let record = bcf::Record::default();
    /// assert_eq!(record.reference_sequence_id()?, 0);
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn reference_sequence_id(&self) -> io::Result<usize> {
        let n = self.0.reference_sequence_id();
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    /// Returns the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::{record::value::{map::Contig, Map}, StringMaps},
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_contig("sq0", Map::<Contig>::new())
    ///     .build();
    /// let string_maps = StringMaps::try_from(&header)?;
    ///
    /// let record = bcf::Record::default();
    /// assert_eq!(record.reference_sequence_name(&string_maps)?, "sq0");
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_sequence_name<'h>(&self, string_maps: &'h StringMaps) -> io::Result<&'h str> {
        self.reference_sequence_id().and_then(|i| {
            string_maps.contigs().get_index(i).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "missing reference sequence name in contig string map",
                )
            })
        })
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
    /// use noodles_core::Position;
    /// let record = bcf::Record::default();
    /// assert_eq!(record.position().transpose()?, Some(Position::MIN));
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn position(&self) -> Option<io::Result<Position>> {
        self.0.position().map(|n| {
            usize::try_from(n)
                .map(|m| m + 1)
                .and_then(Position::try_from)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
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
    /// use noodles_bcf as bcf;
    /// use noodles_core::Position;
    /// let record = bcf::Record::default();
    /// assert_eq!(record.end()?, Position::MIN);
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn end(&self) -> io::Result<Position> {
        let Some(start) = self.position().transpose()? else {
            todo!();
        };

        let len = self.rlen()?;

        start.checked_add(len - 1).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "calculation of the end position overflowed",
            )
        })
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
    /// use noodles_vcf::variant::record::Ids;
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
    /// use noodles_vcf::variant::record::AlternateBases;
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
    /// use noodles_vcf::variant::record::Filters;
    /// let record = bcf::Record::default();
    /// assert!(record.filters().is_empty());
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
    /// use noodles_vcf::variant::record::Info;
    /// let record = bcf::Record::default();
    /// assert!(record.info().is_empty());
    /// ```
    pub fn info(&self) -> Info<'_> {
        self.0.info()
    }

    /// Returns the samples.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// use noodles_vcf::variant::record::Samples;
    /// let record = bcf::Record::default();
    /// assert!(record.samples()?.is_empty());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn samples(&self) -> io::Result<Samples<'_>> {
        self.0.samples()
    }
}

impl fmt::Debug for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Record")
            .field("reference_sequence_id", &self.reference_sequence_id())
            .field("position", &self.position())
            .field("ids", &self.ids())
            .field("reference_bases", &self.reference_bases())
            .field("alternate_bases", &self.alternate_bases())
            .field("quality_score", &self.quality_score())
            .field("filters", &self.filters())
            .field("info", &self.info())
            .field("samples", &self.samples())
            .finish()
    }
}

impl vcf::variant::Record for Record {
    fn reference_sequence_name<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
    ) -> io::Result<&'a str> {
        self.reference_sequence_name(header.string_maps())
    }

    fn position(&self) -> Option<io::Result<Position>> {
        self.position()
    }

    fn ids(&self) -> Box<dyn vcf::variant::record::Ids + '_> {
        Box::new(self.ids())
    }

    fn reference_bases(&self) -> Box<dyn vcf::variant::record::ReferenceBases + '_> {
        Box::new(self.reference_bases())
    }

    fn alternate_bases(&self) -> Box<dyn vcf::variant::record::AlternateBases + '_> {
        Box::new(self.alternate_bases())
    }

    fn quality_score(&self) -> Option<io::Result<f32>> {
        self.quality_score().transpose()
    }

    fn filters(&self) -> Box<dyn vcf::variant::record::Filters + '_> {
        Box::new(self.filters())
    }

    fn info(&self) -> Box<dyn vcf::variant::record::Info + '_> {
        Box::new(self.info())
    }

    fn samples(&self) -> io::Result<Box<dyn vcf::variant::record::Samples + '_>> {
        self.samples()
            .map(|samples| Box::new(samples) as Box<dyn vcf::variant::record::Samples>)
    }
}
