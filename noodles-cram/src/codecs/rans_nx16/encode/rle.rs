mod context;

use std::io;

pub use self::context::{Context, build_context, write_context};
use crate::io::writer::num::write_uint7;

pub fn encode(src: &[u8], ctx: &mut Context) -> io::Result<Vec<u8>> {
    let mut buf = vec![0; src.len()];
    let mut end = 0;

    let mut i = 0;

    while i < src.len() {
        buf[end] = src[i];
        end += 1;

        if ctx.alphabet[src[i] as usize] {
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

    Ok(buf)
}
