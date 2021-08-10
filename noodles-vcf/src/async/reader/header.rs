use std::{
    future::Future,
    pin::Pin,
    str,
    task::{Context, Poll},
};

use futures::ready;
use pin_project_lite::pin_project;
use tokio::io::{self, AsyncBufRead};

const LINE_FEED: u8 = b'\n';

const HEADER_PREFIX: u8 = b'#';

pin_project! {
    pub struct ReadHeader<'a, R> {
        #[pin]
        reader: &'a mut R,
        buf: Vec<u8>,
        i: usize,
        is_eol: bool,
    }
}

impl<'a, R> ReadHeader<'a, R>
where
    R: AsyncBufRead + Unpin,
{
    pub(crate) fn new(reader: &'a mut R) -> Self {
        Self {
            reader,
            buf: Vec::new(),
            i: 0,
            is_eol: false,
        }
    }
}

impl<R> Future for ReadHeader<'_, R>
where
    R: AsyncBufRead + Unpin,
{
    type Output = io::Result<String>;

    fn poll(mut self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Self::Output> {
        let mut this = self.as_mut().project();

        for i in (*this.i).. {
            let buf = ready!(this.reader.as_mut().poll_fill_buf(cx))?;

            if (i == 0 || *this.is_eol) && buf.first().map(|&b| b != HEADER_PREFIX).unwrap_or(true)
            {
                break;
            }

            let (read_eol, len) = if let Some(i) = buf.iter().position(|&b| b == LINE_FEED) {
                this.buf.extend(&buf[..=i]);
                (true, i + 1)
            } else {
                this.buf.extend(buf);
                (false, buf.len())
            };

            *this.i = i;
            *this.is_eol = read_eol;

            this.reader.as_mut().consume(len);
        }

        let result = str::from_utf8(this.buf)
            .map(|s| s.into())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e));

        Poll::Ready(result)
    }
}
