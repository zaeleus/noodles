use std::{
    future::Future,
    pin::Pin,
    task::{Context, Poll},
};

use futures::ready;
use tokio::io::{self, AsyncBufRead};

const LINE_FEED: u8 = b'\n';

const HEADER_PREFIX: u8 = b'#';

pub struct ReadHeader<'a, R> {
    reader: &'a mut R,
}

impl<'a, R> ReadHeader<'a, R>
where
    R: AsyncBufRead + Unpin,
{
    pub(crate) fn new(reader: &'a mut R) -> Self {
        Self { reader }
    }
}

impl<R> Future for ReadHeader<'_, R>
where
    R: AsyncBufRead + Unpin,
{
    type Output = io::Result<String>;

    fn poll(mut self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Self::Output> {
        let reader = Pin::new(&mut self.reader);
        read_header(reader, cx)
    }
}

fn read_header<R>(mut reader: Pin<&mut R>, cx: &mut Context<'_>) -> Poll<io::Result<String>>
where
    R: AsyncBufRead,
{
    let mut header_buf = Vec::new();
    let mut is_eol = false;

    for i in 0.. {
        let buf = ready!(reader.as_mut().poll_fill_buf(cx))?;

        if (i == 0 || is_eol) && buf.first().map(|&b| b != HEADER_PREFIX).unwrap_or(true) {
            break;
        }

        let (read_eol, len) = if let Some(i) = buf.iter().position(|&b| b == LINE_FEED) {
            header_buf.extend(&buf[..=i]);
            (true, i + 1)
        } else {
            header_buf.extend(buf);
            (false, buf.len())
        };

        is_eol = read_eol;

        reader.as_mut().consume(len);
    }

    Poll::Ready(
        String::from_utf8(header_buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    )
}
