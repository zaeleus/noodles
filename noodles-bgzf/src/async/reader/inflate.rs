use std::{
    future::Future,
    io,
    pin::Pin,
    task::{Context, Poll},
};

use bytes::Bytes;
use pin_project_lite::pin_project;
use tokio::task::JoinHandle;

use crate::Block;

pin_project! {
    pub struct Inflate {
        #[pin]
        handle: JoinHandle<io::Result<Block>>,
    }
}

impl Inflate {
    pub(super) fn new(buf: Bytes) -> Self {
        Self {
            handle: tokio::task::spawn_blocking(move || inflate(buf)),
        }
    }
}

impl Future for Inflate {
    type Output = io::Result<Block>;

    fn poll(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Self::Output> {
        self.project().handle.poll(cx)?
    }
}

fn inflate(src: Bytes) -> io::Result<Block> {
    use crate::reader::frame::parse_block;

    let mut block = Block::default();
    parse_block(&src, &mut block)?;
    Ok(block)
}
