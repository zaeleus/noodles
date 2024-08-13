//! Async alignment reader.
//! todo add an example

mod builder;

use futures::{Stream, StreamExt};
use noodles_bam as bam;
use noodles_cram as cram;
use noodles_sam as sam;
use std::pin::Pin;
use tokio::io::{self, AsyncBufRead};

pub use self::builder::Builder;

/// An async alignment reader.
pub enum Reader<R> {
    /// SAM.
    Sam(sam::r#async::io::Reader<R>),
    /// BAM.
    Bam(bam::r#async::io::Reader<R>),
    /// CRAM.
    Cram(cram::r#async::io::Reader<R>),
}

impl<R> Reader<R>
where
    R: AsyncBufRead + Unpin,
{
    /// Reads the SAM header
    ///
    /// # Examples
    /// todo once Builder is implemented
    pub async fn read_header(&mut self) -> io::Result<sam::Header> {
        match self {
            Self::Sam(reader) => reader.read_header().await,
            Self::Bam(reader) => reader.read_header().await,
            Self::Cram(reader) => reader.read_header().await,
        }
    }

    /// Returns an iterator over records starting from the current stream position.
    /// todo add an example
    pub fn records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
    ) -> impl Stream<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'r {
        #[allow(clippy::type_complexity)]
        let records: Pin<
            Box<dyn Stream<Item = io::Result<Box<dyn sam::alignment::Record>>>>,
        > = match self {
            Self::Sam(reader) => Box::pin(
                reader
                    .records()
                    .map(|result| result.map(|r| Box::new(r) as Box<dyn sam::alignment::Record>)),
            ),
            Self::Bam(reader) => Box::pin(
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
