//! Indexed alignment reader.

mod builder;

pub use self::builder::Builder;

use std::io::{self, Read, Seek};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_cram as cram;
use noodles_sam::{self as sam, alignment::RecordBuf};

/// An indexed alignment reader.
pub enum IndexedReader<R> {
    /// SAM.
    Sam(sam::IndexedReader<R>),
    /// BAM.
    Bam(bam::IndexedReader<bgzf::Reader<R>>),
    /// CRAM.
    Cram(cram::IndexedReader<R>),
}

impl<R> IndexedReader<R>
where
    R: Read,
{
    /// Reads the SAM header.
    pub fn read_header(&mut self) -> io::Result<sam::Header> {
        match self {
            Self::Sam(reader) => reader.read_header(),
            Self::Bam(reader) => reader.read_header(),
            Self::Cram(reader) => reader.read_header(),
        }
    }

    /// Returns an iterator over records starting from the current stream position.
    pub fn records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
    ) -> impl Iterator<Item = io::Result<RecordBuf>> + 'r {
        let records: Box<dyn Iterator<Item = io::Result<RecordBuf>>> =
            match self {
                Self::Sam(reader) => Box::new(reader.records(header)),
                Self::Bam(reader) => Box::new(reader.records(header)),
                Self::Cram(reader) => Box::new(reader.records(header).map(|result| {
                    result.and_then(|record| record.try_into_alignment_record(header))
                })),
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
        header: &'h sam::Header,
        region: &Region,
    ) -> io::Result<impl Iterator<Item = io::Result<RecordBuf>> + 'r> {
        let records: Box<dyn Iterator<Item = io::Result<RecordBuf>>> =
            match self {
                Self::Sam(reader) => reader.query(header, region).map(Box::new)?,
                Self::Bam(reader) => reader.query(header, region).map(Box::new)?,
                Self::Cram(reader) => Box::new(reader.query(header, region)?.map(|result| {
                    result.and_then(|record| record.try_into_alignment_record(header))
                })),
            };

        Ok(records)
    }
}
