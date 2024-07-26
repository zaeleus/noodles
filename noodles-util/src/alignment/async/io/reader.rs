//! Async alignment reader.
//! todo add an example

use futures::Stream;
use noodles_bam as bam;
use noodles_cram as cram;
use noodles_sam as sam;
use std::io;
use tokio::io::AsyncBufRead;

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
    /// todo example usage
    pub async fn read_header(&mut self) -> tokio::io::Result<sam::Header> {
        todo!()
    }

    /// Returns an iterator over records starting from the current stream position.
    /// todo add an example
    pub fn records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
    ) -> impl Stream<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'r {
        todo!()
    }
}
