use std::{convert::TryFrom, io};

use noodles_vcf::{self as vcf, record::Position};
use vcf::record::AlternateBases;

use crate::{header::StringMap, reader::record::read_site};

use super::Record;

impl Record {
    pub fn try_into_vcf_record(
        &self,
        header: &vcf::Header,
        string_map: &StringMap,
    ) -> io::Result<vcf::Record> {
        let mut reader = &self[..];
        let (site, _) = read_site(&mut reader, header, string_map)?;

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

        let quality_score = site.qual;
        let filters = site.filter;
        let info = site.info;

        let builder = vcf::Record::builder()
            .set_chromosome(chromosome)
            .set_position(position)
            .set_ids(ids)
            .set_reference_bases(reference_bases)
            .set_alternate_bases(alternate_bases)
            .set_quality_score(quality_score)
            .set_filters(filters)
            .set_info(info);

        builder
            .build()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    }
}
