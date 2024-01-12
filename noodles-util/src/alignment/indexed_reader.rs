//! Indexed alignment reader.

mod builder;

pub use self::builder::Builder;

use std::io::{self, Read, Seek};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_cram as cram;
use noodles_sam::{self as sam, alignment::Record};

/// An indexed alignment reader.
pub enum IndexedReader<R> {
    /// SAM.
    Sam(sam::IndexedReader<R>),
    /// BAM.
    Bam(bam::io::IndexedReader<bgzf::Reader<R>>),
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
    ) -> impl Iterator<Item = io::Result<Box<dyn Record>>> + 'r {
        let records: Box<dyn Iterator<Item = io::Result<Box<dyn Record>>>> = match self {
            Self::Sam(reader) => Box::new(
                reader
                    .record_bufs(header)
                    .map(|result| result.map(|record| Box::new(record) as Box<dyn Record>)),
            ),
            Self::Bam(reader) => Box::new(
                reader
                    .record_bufs(header)
                    .map(|result| result.map(|record| Box::new(record) as Box<dyn Record>)),
            ),
            Self::Cram(reader) => Box::new(reader.records(header).map(|result| {
                result.and_then(|record| {
                    record
                        .try_into_alignment_record(header)
                        .map(|alignment_record| {
                            Box::new(alignment_record) as Box<dyn sam::alignment::Record>
                        })
                })
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
    ) -> io::Result<impl Iterator<Item = io::Result<Box<dyn Record>>> + 'r> {
        let records: Box<dyn Iterator<Item = io::Result<Box<dyn Record>>>> = match self {
            Self::Sam(reader) => {
                let query = reader.query(header, region)?;

                Box::new(
                    query.map(|result| result.map(|record| Box::new(record) as Box<dyn Record>)),
                )
            }
            Self::Bam(reader) => {
                let query = reader.query(header, region)?;

                Box::new(
                    query.map(|result| result.map(|record| Box::new(record) as Box<dyn Record>)),
                )
            }
            Self::Cram(reader) => {
                let query = reader.query(header, region)?;

                Box::new(query.map(|result| {
                    result.and_then(|record| {
                        record
                            .try_into_alignment_record(header)
                            .map(|alignment_record| {
                                Box::new(alignment_record) as Box<dyn sam::alignment::Record>
                            })
                    })
                }))
            }
        };

        Ok(records)
    }
}
