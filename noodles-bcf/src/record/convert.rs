use std::io;

use noodles_vcf::{self as vcf, record::Position};
use vcf::record::{AlternateBases, Format, QualityScore};

use crate::{header::StringMap, reader::record::read_record};

use super::{value::Float, Record};

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
    /// let record = bcf::Record::from(vec![
    ///     0x00, 0x00, 0x00, 0x00, // chrom = sq0
    ///     0x00, 0x00, 0x00, 0x00, // pos = 0 (base 0)
    ///     0x01, 0x00, 0x00, 0x00, // rlen = 1
    ///     0x01, 0x00, 0x80, 0x7f, // qual = Float::Missing
    ///     0x00, 0x00, // n_info = 0
    ///     0x01, 0x00, // n_allele = 1
    ///     0x00, // n_sample = 0
    ///     0x00, 0x00, 0x00, // n_fmt = 0
    ///     0x07, // id = [missing]
    ///     0x17, 0x41, // ref = A
    ///     0x00, // filter = []
    /// ]);
    ///
    /// let actual = record.try_into_vcf_record(&header, &string_map)?;
    ///
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
        let mut reader = &self[..];
        let (site, genotypes) = read_record(&mut reader, header, string_map)?;

        let (_, contig) = usize::try_from(site.chrom)
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

        let position = Position::try_from(site.pos + 1)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let quality_score = match site.qual {
            Float::Value(value) => QualityScore::try_from(value)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?,
            Float::Missing => QualityScore::default(),
            qual => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid qual: {:?}", qual),
                ));
            }
        };

        let ids = site.id;

        let (raw_reference_bases, raw_alternate_bases) = site.ref_alt.split_at(1);

        let reference_bases = raw_reference_bases
            .first()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing reference bases"))
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            })?;

        let alternate_alleles: Vec<_> = raw_alternate_bases
            .iter()
            .map(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            })
            .collect::<Result<_, _>>()?;

        let alternate_bases = AlternateBases::from(alternate_alleles);

        let info = site.info;

        let mut builder = vcf::Record::builder()
            .set_chromosome(chromosome)
            .set_position(position)
            .set_ids(ids)
            .set_reference_bases(reference_bases)
            .set_alternate_bases(alternate_bases)
            .set_quality_score(quality_score)
            .set_info(info);

        if let Some(filters) = site.filter {
            builder = builder.set_filters(filters);
        }

        if let Some(first_genotype) = genotypes.first() {
            let keys: Vec<_> = first_genotype.keys().cloned().collect();
            let format = Format::try_from(keys)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
            builder = builder.set_format(format).set_genotypes(genotypes);
        }

        builder
            .build()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    }
}
