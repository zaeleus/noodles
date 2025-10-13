use std::{io, num::NonZero};

use crate::io::reader::num::{read_u8, read_uint7_as};

pub struct Context<'a> {
    pub symbol_count: NonZero<usize>,
    pub mapping_table: &'a [u8],
    pub uncompressed_size: usize,
}

pub fn read_context<'a>(
    src: &mut &'a [u8],
    uncompressed_size: usize,
) -> io::Result<(Context<'a>, usize)> {
    let symbol_count = read_u8(src).and_then(|n| {
        NonZero::try_from(usize::from(n)).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let (mapping_table, rest) = src
        .split_at_checked(symbol_count.get())
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    let len = read_uint7_as(src)?;

    let context = Context {
        symbol_count,
        mapping_table,
        uncompressed_size,
    };

    Ok((context, len))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_context() -> io::Result<()> {
        let src = [
            0x04, // symbol count = 4
            b'A', b'C', b'G', b'T', // mapping table
            0x08, // uncompressed size = 8
        ];

        let (ctx, len) = read_context(&mut &src[..], 13)?;

        assert_eq!(ctx.symbol_count, const { NonZero::new(4).unwrap() });
        assert_eq!(ctx.mapping_table, b"ACGT");
        assert_eq!(ctx.uncompressed_size, 13);
        assert_eq!(len, 8);

        Ok(())
    }
}
