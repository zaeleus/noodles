use std::pin::Pin;

use futures::{Stream, StreamExt};
use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf as vcf;
use tokio::io::{self, AsyncRead, BufReader};

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
