use std::io;

use noodles_vcf as vcf;

use super::Record;
use crate::header::StringMaps;

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
    /// let string_maps = raw_header.parse()?;
    ///
    /// let record = bcf::lazy::Record::default();
    ///
    /// let actual = record.try_into_vcf_record(&header, &string_maps)?;
    /// let expected = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("N")
    ///     .build()?;
    ///
    /// assert_eq!(actual, expected);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_into_vcf_record(
        &self,
        header: &vcf::Header,
        string_maps: &StringMaps,
    ) -> io::Result<vcf::Record> {
        let chromosome = string_maps
            .contigs()
            .get_index(self.chromosome_id())
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid chrom"))?;

        let filters = self
            .filters()
            .try_into_vcf_record_filters(string_maps.strings())?;

        let info = self
            .info()
            .try_into_vcf_record_info(header, string_maps.strings())?;

        let genotypes = self
            .genotypes()
            .try_into_vcf_record_genotypes(header, string_maps.strings())?;

        let mut builder = vcf::Record::builder()
            .set_chromosome(chromosome)
            .set_position(self.position())
            .set_ids(self.ids().clone())
            .set_reference_bases(self.reference_bases())
            .set_alternate_bases(self.alternate_bases().clone())
            .set_info(info)
            .set_genotypes(genotypes);

        if let Some(quality_score) = self.quality_score() {
            builder = builder.set_quality_score(quality_score);
        }

        if let Some(filters) = filters {
            builder = builder.set_filters(filters);
        }

        builder
            .build()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    }
}
