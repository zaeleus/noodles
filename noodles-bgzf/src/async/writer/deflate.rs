use std::{
    future::Future,
    io::{self, Write},
    pin::Pin,
    task::{Context, Poll},
};

use bytes::BytesMut;
use flate2::{write::DeflateEncoder, Compression, Crc};
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
    let mut encoder = DeflateEncoder::new(Vec::new(), compression);
    encoder.write_all(&data[..])?;
    let compressed_data = encoder.finish()?;

    let mut crc = Crc::new();
    crc.update(&data[..]);

    Ok((compressed_data, crc.sum(), crc.amount()))
}
