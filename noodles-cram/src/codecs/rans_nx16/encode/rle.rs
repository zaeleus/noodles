mod context;

use std::io;

pub use self::context::{Context, build_context, write_context};
use crate::io::writer::num::write_uint7;

pub fn encode(src: &[u8], ctx: &mut Context) -> io::Result<Vec<u8>> {
    let mut iter = src.iter().peekable();
    let mut dst = Vec::new();

    while let Some(&sym) = iter.next() {
        write_u8(&mut dst, sym);

        if ctx.alphabet[usize::from(sym)] {
            let mut len = 0;

            while let Some(&&s) = iter.peek()
                && s == sym
            {
                len += 1;
                iter.next();
            }

            let n =
                u32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
            write_uint7(&mut ctx.dst, n)?;
        }
    }

    Ok(dst)
}

fn write_u8(dst: &mut Vec<u8>, n: u8) {
    dst.push(n);
}
