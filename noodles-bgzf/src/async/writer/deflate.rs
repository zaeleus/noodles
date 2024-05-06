use std::{
    future::Future,
    io,
    pin::Pin,
    task::{Context, Poll},
};

use bytes::BytesMut;
use pin_project_lite::pin_project;
use tokio::task::JoinHandle;

use super::CompressionLevel;
use crate::deflate;

// (CDATA, CRC32, ISIZE)
pub type GzData = (Vec<u8>, u32, usize);

pin_project! {
    pub struct Deflate {
        #[pin]
        handle: JoinHandle<io::Result<GzData>>,
    }
}

impl Deflate {
    pub fn new(data: BytesMut, compression_level: CompressionLevel) -> Self {
        Self {
            handle: tokio::task::spawn_blocking(move || {
                let mut dst = Vec::new();

                deflate::encode(&data, compression_level, &mut dst)
                    .map(|crc32| (dst, crc32, data.len()))
            }),
        }
    }
}

impl Future for Deflate {
    type Output = io::Result<GzData>;

    fn poll(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Self::Output> {
        self.project().handle.poll(cx)?
    }
}
