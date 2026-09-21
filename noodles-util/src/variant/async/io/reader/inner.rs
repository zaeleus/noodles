use std::pin::Pin;

use futures::{Stream, StreamExt};
use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_vcf as vcf;
use tokio::io::{self, AsyncRead, AsyncSeek, BufReader};

use crate::variant::{Index, Record};

pub(super) enum Inner<R>
where
    R: AsyncRead,
{
    Bcf(bcf::r#async::io::Reader<bgzf::r#async::io::Reader<BufReader<R>>>),
    BcfRaw(bcf::r#async::io::Reader<BufReader<R>>),
    Vcf(vcf::r#async::io::Reader<BufReader<R>>),
    VcfGz(vcf::r#async::io::Reader<bgzf::r#async::io::Reader<BufReader<R>>>),
}

impl<R> Inner<R>
where
    R: AsyncRead + Unpin,
{
    pub(super) async fn read_header(&mut self) -> io::Result<vcf::Header> {
        match self {
            Self::Bcf(reader) => reader.read_header().await,
            Self::BcfRaw(reader) => reader.read_header().await,
            Self::Vcf(reader) => reader.read_header().await,
            Self::VcfGz(reader) => reader.read_header().await,
        }
    }

    pub(super) async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        match self {
            Self::Bcf(reader) => {
                if !matches!(record, Record::Bcf(_)) {
                    *record = Record::Bcf(bcf::Record::default());
                }

                if let Record::Bcf(r) = record {
                    reader.read_record(r).await
                } else {
                    unreachable!();
                }
            }
            Self::BcfRaw(reader) => {
                if !matches!(record, Record::Bcf(_)) {
                    *record = Record::Bcf(bcf::Record::default());
                }

                if let Record::Bcf(r) = record {
                    reader.read_record(r).await
                } else {
                    unreachable!();
                }
            }
            Self::Vcf(reader) => {
                if !matches!(record, Record::Vcf(_)) {
                    *record = Record::Vcf(vcf::Record::default());
                }

                if let Record::Vcf(r) = record {
                    reader.read_record(r).await
                } else {
                    unreachable!();
                }
            }
            Self::VcfGz(reader) => {
                if !matches!(record, Record::Vcf(_)) {
                    *record = Record::Vcf(vcf::Record::default());
                }

                if let Record::Vcf(r) = record {
                    reader.read_record(r).await
                } else {
                    unreachable!();
                }
            }
        }
    }

    pub(super) fn records(
        &mut self,
    ) -> impl Stream<Item = io::Result<Box<dyn vcf::variant::Record>>> + '_ {
        let records: Pin<Box<dyn Stream<Item = io::Result<_>>>> = match self {
            Self::Bcf(reader) => Box::pin(reader.records().map(|result| {
                result.map(|record| Box::new(record) as Box<dyn vcf::variant::Record>)
            })),
            Self::BcfRaw(reader) => Box::pin(reader.records().map(|result| {
                result.map(|record| Box::new(record) as Box<dyn vcf::variant::Record>)
            })),
            Self::Vcf(reader) => Box::pin(reader.records().map(|result| {
                result.map(|record| Box::new(record) as Box<dyn vcf::variant::Record>)
            })),
            Self::VcfGz(reader) => Box::pin(reader.records().map(|result| {
                result.map(|record| Box::new(record) as Box<dyn vcf::variant::Record>)
            })),
        };

        records
    }
}

impl<R> Inner<R>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    pub(super) fn query<'r, 'h: 'r>(
        &'r mut self,
        header: &'h vcf::Header,
        index: &Index,
        region: &Region,
    ) -> io::Result<impl Stream<Item = io::Result<Box<dyn vcf::variant::Record>>> + 'r> {
        let records: Pin<Box<dyn Stream<Item = io::Result<_>>>> = match (self, index) {
            (Self::VcfGz(reader), Index::Vcf(idx)) => {
                let query = reader.query(header, idx, region)?;

                Box::pin(query.records().map(|result| {
                    result.map(|record| Box::new(record) as Box<dyn vcf::variant::Record>)
                }))
            }
            (Self::Bcf(reader), Index::Bcf(idx)) => {
                let query = reader.query(header, idx, region)?;

                Box::pin(query.records().map(|result| {
                    result.map(|record| Box::new(record) as Box<dyn vcf::variant::Record>)
                }))
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "format-index mismatch",
                ));
            }
        };

        Ok(records)
    }
}
