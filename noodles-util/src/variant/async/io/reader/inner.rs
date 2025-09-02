use std::pin::Pin;

use futures::{Stream, StreamExt};
use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf as vcf;
use tokio::io::{self, AsyncRead, BufReader};

use crate::variant::Record;

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
            Inner::Bcf(reader) => {
                if let Record::Bcf(bcf) = record {
                    reader.read_record(bcf).await
                } else {
                    let mut bcf = bcf::Record::default();
                    let ret = reader.read_record(&mut bcf).await?;
                    *record = Record::Bcf(bcf);
                    Ok(ret)
                }
            }
            Inner::BcfRaw(reader) => {
                if let Record::Bcf(bcf) = record {
                    reader.read_record(bcf).await
                } else {
                    let mut bcf = bcf::Record::default();
                    let ret = reader.read_record(&mut bcf).await?;
                    *record = Record::Bcf(bcf);
                    Ok(ret)
                }
            }
            Inner::Vcf(reader) => {
                if let Record::Vcf(vcf) = record {
                    reader.read_record(vcf).await
                } else {
                    let mut vcf = vcf::Record::default();
                    let ret = reader.read_record(&mut vcf).await?;
                    *record = Record::Vcf(vcf);
                    Ok(ret)
                }
            }
            Inner::VcfGz(reader) => {
                if let Record::Vcf(vcf) = record {
                    reader.read_record(vcf).await
                } else {
                    let mut vcf = vcf::Record::default();
                    let ret = reader.read_record(&mut vcf).await?;
                    *record = Record::Vcf(vcf);
                    Ok(ret)
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
