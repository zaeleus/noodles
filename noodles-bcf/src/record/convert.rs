use std::{io, str};

use noodles_vcf::{
    self as vcf,
    header::StringMaps,
    variant::record::{AlternateBases, Ids, Info},
};

use super::Record;

impl Record {
    /// Converts a VCF record to a BCF record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// use noodles_core::Position;
    /// use noodles_vcf as vcf;
    ///
    /// let raw_header = "##fileformat=VCFv4.3\n##contig=<ID=sq0>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    /// let header: vcf::Header = raw_header.parse()?;
    /// let string_maps = raw_header.parse()?;
    ///
    /// let record = bcf::Record::default();
    ///
    /// let actual = record.try_into_vcf_record(&header, &string_maps)?;
    /// let expected = vcf::variant::RecordBuf::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_position(Position::MIN)
    ///     .set_reference_bases("N")
    ///     .build();
    ///
    /// assert_eq!(actual, expected);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_into_vcf_record(
        &self,
        header: &vcf::Header,
        string_maps: &StringMaps,
    ) -> io::Result<vcf::variant::RecordBuf> {
        let chromosome = string_maps
            .contigs()
            .get_index(self.reference_sequence_id()?)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid chrom"))?;

        let raw_reference_bases = self.reference_bases();
        let reference_bases = str::from_utf8(raw_reference_bases.as_ref())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let alternate_bases: Vec<_> = self
            .alternate_bases()
            .iter()
            .map(|result| result.map(String::from))
            .collect::<io::Result<_>>()?;

        let info = self
            .info()
            .iter(header)
            .map(|result| {
                result.and_then(|(key, value)| {
                    let v = value.map(|v| v.try_into()).transpose()?;
                    Ok((key.into(), v))
                })
            })
            .collect::<io::Result<_>>()?;

        let samples = self.samples()?.try_into_vcf_record_samples(header)?;

        let mut builder = vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name(chromosome)
            .set_reference_bases(reference_bases)
            .set_alternate_bases(alternate_bases.into())
            .set_info(info)
            .set_samples(samples);

        if let Some(position) = self.position().transpose()? {
            builder = builder.set_position(position);
        }

        if !self.ids().is_empty() {
            const DELIMITER: char = ';';

            let ids = str::from_utf8(self.ids().as_ref())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
                .map(|s| s.split(DELIMITER).map(String::from).collect())?;

            builder = builder.set_ids(ids);
        }

        if let Some(quality_score) = self.quality_score()? {
            builder = builder.set_quality_score(quality_score);
        }

        if !self.filters().is_empty()? {
            let filters = self
                .filters()
                .iter(header)?
                .map(|filter| filter.map(String::from))
                .collect::<io::Result<_>>()?;

            builder = builder.set_filters(filters);
        }

        Ok(builder.build())
    }
}
