use std::{
    future::Future,
    pin::Pin,
    task::{Context, Poll},
};

use futures::Stream;
use tokio::{
    io::{self, AsyncRead},
    pin,
};

use crate::Record;

use super::Reader;

pub struct Records<'a, R>
where
    R: AsyncRead,
{
    reader: &'a mut Reader<R>,
    record: Record,
}

impl<'a, R> Records<'a, R>
where
    R: AsyncRead + Unpin,
{
    pub(crate) fn new(reader: &'a mut Reader<R>) -> Self {
        Self {
            reader,
            record: Record::default(),
        }
    }

    async fn read_record(&mut self) -> Option<io::Result<Record>> {
        match self.reader.read_record(&mut self.record).await {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}

impl<'a, R> Stream for Records<'a, R>
where
    R: AsyncRead + Unpin,
{
    type Item = io::Result<Record>;

    fn poll_next(mut self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        let read_record = self.read_record();
        pin!(read_record);
        read_record.poll(cx)
    }
}
