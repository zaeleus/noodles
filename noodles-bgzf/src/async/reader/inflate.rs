use std::{
    future::Future,
    io,
    pin::Pin,
    task::{Context, Poll},
};

use bytes::{Buf, Bytes};
use pin_project_lite::pin_project;
use tokio::task::JoinHandle;

use crate::{gz, Block};

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
    if !is_valid_header(&mut src) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BGZF header",
        ));
    }

    let bsize = u64::from(src.get_u16_le()) + 1;

    let cdata = src.split_to(src.len() - gz::TRAILER_SIZE);

    // trailer
    let crc32 = src.get_u32_le();
    let r#isize = usize::try_from(src.get_u32_le())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let mut block = Block::default();

    block.set_size(bsize);

    let data = block.data_mut();
    data.set_position(0);
    data.resize(r#isize);

    crate::reader::block::inflate(&cdata, crc32, data.as_mut())?;

    Ok(block)
}

fn is_valid_header(src: &mut Bytes) -> bool {
    use std::mem;

    const BGZF_CM: u8 = 0x08; // DEFLATE
    const BGZF_FLG: u8 = 0x04; // FEXTRA
    const BGZF_XLEN: u16 = 6;
    const BGZF_SI1: u8 = b'B';
    const BGZF_SI2: u8 = b'C';
    const BGZF_SLEN: u16 = 2;

    let magic_number = src.split_to(gz::MAGIC_NUMBER.len());
    let cm = src.get_u8();
    let flg = src.get_u8();

    // 4 (MTIME) + 1 (XFL) + 1 (OS)
    src.advance(mem::size_of::<u32>() + mem::size_of::<u8>() + mem::size_of::<u8>());

    let xlen = src.get_u16_le();
    let subfield_id_1 = src.get_u8();
    let subfield_id_2 = src.get_u8();
    let bc_len = src.get_u16_le();

    magic_number[..] == gz::MAGIC_NUMBER
        && cm == BGZF_CM
        && flg == BGZF_FLG
        && xlen == BGZF_XLEN
        && subfield_id_1 == BGZF_SI1
        && subfield_id_2 == BGZF_SI2
        && bc_len == BGZF_SLEN
}
