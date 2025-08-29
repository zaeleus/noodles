use std::pin::Pin;

use futures::{Stream, StreamExt};
use noodles_bcf as bcf;
use noodles_vcf as vcf;
use tokio::io::{self, AsyncBufRead};

pub(super) enum Inner<R> {
    Vcf(vcf::r#async::io::Reader<R>),
    Bcf(bcf::r#async::io::Reader<R>),
}

impl<R> Inner<R>
where
    R: AsyncBufRead + Unpin,
{
    pub(super) async fn read_header(&mut self) -> io::Result<vcf::Header> {
        match self {
            Self::Bcf(reader) => reader.read_header().await,
            Self::Vcf(reader) => reader.read_header().await,
        }
    }

    pub(super) fn records(
        &mut self,
    ) -> impl Stream<Item = io::Result<Box<dyn vcf::variant::Record>>> + '_ {
        #[allow(clippy::type_complexity)]
        let records: Pin<
            Box<dyn Stream<Item = io::Result<Box<dyn vcf::variant::Record>>>>,
        > = match self {
            Self::Bcf(reader) => Box::pin(reader.records().map(|result| {
                result.map(|record| Box::new(record) as Box<dyn vcf::variant::Record>)
            })),
            Self::Vcf(reader) => Box::pin(reader.records().map(|result| {
                result.map(|record| Box::new(record) as Box<dyn vcf::variant::Record>)
            })),
        };

        records
    }
}
