use std::{convert::TryFrom, io};

use noodles_vcf::{self as vcf, record::Position};
use vcf::record::{AlternateBases, Format, QualityScore};

use crate::{header::StringMap, reader::record::read_record};

use super::{value::Float, Record};

impl Record {
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
            qual => todo!("unhandled quality score value: {:?}", qual),
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

        let filters = site.filter;
        let info = site.info;

        let mut builder = vcf::Record::builder()
            .set_chromosome(chromosome)
            .set_position(position)
            .set_ids(ids)
            .set_reference_bases(reference_bases)
            .set_alternate_bases(alternate_bases)
            .set_quality_score(quality_score)
            .set_filters(filters)
            .set_info(info);

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
