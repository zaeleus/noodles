//! Indexed variant reader.

mod builder;

pub use self::builder::Builder;

use std::io::{self, BufRead};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_vcf::{self as vcf, variant::Record};

/// An indexed variant reader.
pub enum IndexedReader<R> {
    /// VCF.
    Vcf(vcf::io::IndexedReader<R>),
    /// BCF.
    Bcf(bcf::io::IndexedReader<R>),
}

impl<R> IndexedReader<R>
where
    R: BufRead,
{
    /// Reads the VCF header.
    pub fn read_header(&mut self) -> io::Result<vcf::Header> {
        match self {
            Self::Vcf(reader) => reader.read_header(),
            Self::Bcf(reader) => reader.read_header(),
        }
    }

    /// Returns an iterator over records starting from the current stream position.
    pub fn records(&mut self) -> impl Iterator<Item = io::Result<Box<dyn Record>>> {
        let records: Box<dyn Iterator<Item = io::Result<Box<dyn Record>>>> = match self {
            Self::Vcf(reader) => Box::new(
                reader
                    .records()
                    .map(|result| result.map(|record| Box::new(record) as Box<dyn Record>)),
            ),
            Self::Bcf(reader) => Box::new(
                reader
                    .records()
                    .map(|result| result.map(|record| Box::new(record) as Box<dyn Record>)),
            ),
        };

        records
    }
}

impl<R> IndexedReader<R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    /// Returns an iterator over records that intersects the given region.
    pub fn query<'r, 'h: 'r>(
        &'r mut self,
        header: &'h vcf::Header,
        region: &Region,
    ) -> io::Result<impl Iterator<Item = io::Result<Box<dyn Record>>> + use<'r, 'h, R>> {
        let records: Box<dyn Iterator<Item = io::Result<Box<dyn Record>>>> = match self {
            Self::Vcf(reader) => Box::new(
                reader
                    .query(header, region)?
                    .map(|result| result.map(|record| Box::new(record) as Box<dyn Record>)),
            ),
            Self::Bcf(reader) => Box::new(
                reader
                    .query(header, region)?
                    .map(|result| result.map(|record| Box::new(record) as Box<dyn Record>)),
            ),
        };

        Ok(records)
    }
}
