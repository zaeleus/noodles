use std::{
    future::Future,
    io,
    pin::Pin,
    task::{Context, Poll},
};

use bytes::{Buf, Bytes};
use flate2::Crc;
use pin_project_lite::pin_project;
use tokio::task::JoinHandle;

use crate::{gz, Block, BGZF_HEADER_SIZE};

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

fn inflate(mut src: Bytes) -> io::Result<Block> {
    use crate::reader::inflate_data;

    let mut header = src.split_to(BGZF_HEADER_SIZE);
    header.advance(16); // [ID1, ..., SLEN]
    let bsize = u64::from(header.get_u16_le()) + 1;

    let cdata = src.split_to(src.len() - gz::TRAILER_SIZE);

    // trailer
    let crc32 = src.get_u32_le();
    let r#isize = usize::try_from(src.get_u32_le())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let mut block = Block::default();

    block.set_clen(bsize);
    block.set_upos(0);
    block.set_ulen(r#isize);

    inflate_data(&cdata, block.buffer_mut())?;

    let mut crc = Crc::new();
    crc.update(block.buffer());

    if crc.sum() == crc32 {
        Ok(block)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "block data checksum mismatch",
        ))
    }
}
