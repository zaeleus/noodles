use std::{convert::TryFrom, io};

use noodles_vcf::{self as vcf, record::Position};

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

        let reference_bases = site
            .ref_alt
            .first()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing reference bases"))
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            })?;

        let builder = vcf::Record::builder()
            .set_chromosome(chromosome)
            .set_position(position)
            .set_reference_bases(reference_bases);

        builder
            .build()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    }
}
