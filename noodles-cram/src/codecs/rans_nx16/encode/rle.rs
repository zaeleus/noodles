mod context;

use std::io::{self, Write};

pub use self::context::{Context, build_context};
use crate::io::writer::num::write_uint7;

pub fn encode(src: &[u8], ctx: &mut Context) -> io::Result<(Vec<u8>, Vec<u8>)> {
    let mut buf = vec![0; src.len()];
    let mut end = 0;

    let mut i = 0;

    while i < src.len() {
        buf[end] = src[i];
        end += 1;

        if ctx.alphabet[src[i] as usize] > 0 {
            let mut run = 0;
            let last = src[i];

            while i + run + 1 < src.len() && src[i + run + 1] == last {
                run += 1;
            }

            write_uint7(&mut ctx.dst, run as u32)?;

            i += run;
        }

        i += 1;
    }

    buf.truncate(end);

    let mut header = Vec::new();

    let rle_meta_len = ctx.dst.len() as u32;
    write_uint7(&mut header, (rle_meta_len << 1) | 1)?;

    let len = buf.len() as u32;
    write_uint7(&mut header, len)?;

    header.write_all(&ctx.dst)?;

    Ok((header, buf))
}
