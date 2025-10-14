use std::{borrow::Cow, io};

use crate::{
    codecs::rans_nx16::decode::{order_0, split_off},
    io::reader::num::{read_uint7, read_uint7_as},
};

pub struct Context<'a> {
    pub src: Cow<'a, [u8]>,
    pub len: usize,
}

pub fn read_context<'a>(
    src: &mut &'a [u8],
    state_count: usize,
    uncompressed_size: usize,
) -> io::Result<(Context<'a>, usize)> {
    let (context_size, is_compressed) = read_header(src)?;

    let len = read_uint7_as(src)?;
    let context_src = read_src(src, state_count, context_size, is_compressed)?;

    let ctx = Context {
        src: context_src,
        len: uncompressed_size,
    };

    Ok((ctx, len))
}

fn read_header(src: &mut &[u8]) -> io::Result<(usize, bool)> {
    let n = read_uint7(src)?;

    let context_size =
        usize::try_from(n >> 1).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let is_compressed = (n & 0x01) == 0;

    Ok((context_size, is_compressed))
}

fn read_src<'a>(
    src: &mut &'a [u8],
    state_count: usize,
    len: usize,
    is_compressed: bool,
) -> io::Result<Cow<'a, [u8]>> {
    if is_compressed {
        let compressed_size = read_uint7_as(src)?;
        let mut buf = split_off(src, compressed_size)?;

        let mut dst = vec![0; len];
        order_0::decode(&mut buf, &mut dst, state_count)?;
        Ok(Cow::from(dst))
    } else {
        split_off(src, len).map(Cow::from)
    }
}
