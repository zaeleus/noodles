//! Indexed variant reader.

mod builder;

pub use self::builder::Builder;

use std::io::{self, Read, Seek};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_vcf::{self as vcf, Record};

/// An indexed variant reader.
pub enum IndexedReader<R> {
    /// VCF.
    Vcf(vcf::io::IndexedReader<R>),
    /// BCF.
    Bcf(bcf::IndexedReader<bgzf::Reader<R>>),
}

impl<R> IndexedReader<R>
where
    R: Read,
{
    /// Reads the VCF header.
    pub fn read_header(&mut self) -> io::Result<vcf::Header> {
        match self {
            Self::Vcf(reader) => reader.read_header(),
            Self::Bcf(reader) => reader.read_header(),
        }
    }

    /// Returns an iterator over records starting from the current stream position.
    pub fn records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h vcf::Header,
    ) -> impl Iterator<Item = io::Result<Record>> + '_ {
        let records: Box<dyn Iterator<Item = io::Result<Record>>> = match self {
            Self::Vcf(reader) => Box::new(reader.records(header)),
            Self::Bcf(reader) => Box::new(reader.records(header)),
        };

        records
    }
}

impl<R> IndexedReader<R>
where
    R: Read + Seek,
{
    /// Returns an iterator over records that intersects the given region.
    pub fn query<'r, 'h: 'r>(
        &'r mut self,
        header: &'h vcf::Header,
        region: &Region,
    ) -> io::Result<impl Iterator<Item = io::Result<Record>> + '_> {
        let records: Box<dyn Iterator<Item = io::Result<Record>>> = match self {
            Self::Vcf(reader) => reader.query(header, region).map(Box::new)?,
            Self::Bcf(reader) => reader.query(header, region).map(Box::new)?,
        };

        Ok(records)
    }
}
