use std::io;

use noodles_vcf as vcf;

use super::Record;
use crate::header::StringMap;

impl Record {
    /// Converts a VCF record to a BCF record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// use noodles_vcf::{self as vcf, record::Position};
    ///
    /// let raw_header = "##fileformat=VCFv4.3\n##contig=<ID=sq0>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    /// let header: vcf::Header = raw_header.parse()?;
    /// let string_map = raw_header.parse()?;
    ///
    /// let record = bcf::Record::default();
    ///
    /// let actual = record.try_into_vcf_record(&header, &string_map)?;
    /// let expected = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(actual, expected);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_into_vcf_record(
        &self,
        header: &vcf::Header,
        string_map: &StringMap,
    ) -> io::Result<vcf::Record> {
        let (_, contig) = usize::try_from(self.chromosome_id())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            .and_then(|i| {
                header
                    .contigs()
                    .get_index(i)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid chrom"))
            })?;

        let chromosome = contig
            .id()
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let filters = self.filters().try_into_vcf_record_filters(string_map)?;
        let info = self.info().try_into_vcf_record_info(header, string_map)?;

        let mut builder = vcf::Record::builder()
            .set_chromosome(chromosome)
            .set_position(self.position())
            .set_ids(self.ids().clone())
            .set_reference_bases(self.reference_bases().clone())
            .set_alternate_bases(self.alternate_bases().clone())
            .set_info(info);

        if let Some(quality_score) = self.quality_score() {
            builder = builder.set_quality_score(quality_score);
        }

        if let Some(filters) = filters {
            builder = builder.set_filters(filters);
        }

        let (format, genotypes) = self
            .genotypes()
            .try_into_vcf_record_genotypes(header, string_map)?;

        if let Some(format) = format {
            builder = builder.set_format(format).set_genotypes(genotypes);
        }

        builder
            .build()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    }
}
