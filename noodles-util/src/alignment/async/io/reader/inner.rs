use std::pin::Pin;

use futures::{Stream, StreamExt};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_sam as sam;
use tokio::io::{self, AsyncRead, BufReader};

pub(super) enum Inner<R>
where
    R: AsyncRead,
{
    Sam(sam::r#async::io::Reader<BufReader<R>>),
    SamGz(sam::r#async::io::Reader<bgzf::r#async::io::Reader<BufReader<R>>>),
    Bam(bam::r#async::io::Reader<bgzf::r#async::io::Reader<BufReader<R>>>),
    BamRaw(bam::r#async::io::Reader<BufReader<R>>),
    Cram(cram::r#async::io::Reader<BufReader<R>>),
}

impl<R> Inner<R>
where
    R: AsyncRead + Unpin,
{
    pub(super) async fn read_header(&mut self) -> io::Result<sam::Header> {
        match self {
            Self::Sam(reader) => reader.read_header().await,
            Self::SamGz(reader) => reader.read_header().await,
            Self::Bam(reader) => reader.read_header().await,
            Self::BamRaw(reader) => reader.read_header().await,
            Self::Cram(reader) => reader.read_header().await,
        }
    }

    pub(super) fn records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
    ) -> impl Stream<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'r {
        #[allow(clippy::type_complexity)]
        let records: Pin<Box<dyn Stream<Item = io::Result<_>>>> = match self {
            Self::Sam(reader) => Box::pin(
                reader
                    .records()
                    .map(|result| result.map(|r| Box::new(r) as Box<dyn sam::alignment::Record>)),
            ),
            Self::SamGz(reader) => Box::pin(
                reader
                    .records()
                    .map(|result| result.map(|r| Box::new(r) as Box<dyn sam::alignment::Record>)),
            ),
            Self::Bam(reader) => Box::pin(
                reader
                    .records()
                    .map(|result| result.map(|r| Box::new(r) as Box<dyn sam::alignment::Record>)),
            ),
            Self::BamRaw(reader) => Box::pin(
                reader
                    .records()
                    .map(|result| result.map(|r| Box::new(r) as Box<dyn sam::alignment::Record>)),
            ),

            Self::Cram(reader) => Box::pin(
                reader
                    .records(header)
                    .map(|result| result.map(|r| Box::new(r) as Box<dyn sam::alignment::Record>)),
            ),
        };

        records
    }
}
