use std::{
    future::Future,
    io,
    pin::Pin,
    task::{Context, Poll},
};

use bytes::BytesMut;
use flate2::Compression;
use pin_project_lite::pin_project;
use tokio::task::JoinHandle;

// (CDATA, CRC32, ISIZE)
pub type GzData = (Vec<u8>, u32, u32);

pin_project! {
    pub struct Deflate {
        #[pin]
        handle: JoinHandle<io::Result<GzData>>,
    }
}

impl Deflate {
    pub fn new(data: BytesMut, compression: Compression) -> Self {
        Self {
            handle: tokio::task::spawn_blocking(move || deflate(data, compression)),
        }
    }
}

impl Future for Deflate {
    type Output = io::Result<GzData>;

    fn poll(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Self::Output> {
        self.project().handle.poll(cx)?
    }
}

fn deflate(data: BytesMut, compression: Compression) -> io::Result<GzData> {
    use crate::writer::deflate_data;
    deflate_data(&data, compression)
}
