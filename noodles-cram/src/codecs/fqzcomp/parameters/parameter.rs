mod flags;

pub use self::flags::Flags;

use std::{io, num::NonZero};

use super::read_array;
use crate::io::reader::num::{read_u8, read_u16_le};

pub struct Parameter {
    pub context: u16,
    pub flags: Flags,
    pub symbol_count: NonZero<usize>,
    pub q_bits: u8,
    pub q_shift: u8,
    pub q_loc: u8,
    pub s_loc: u8,
    pub p_loc: u8,
    pub d_loc: u8,
    pub q_map: Option<Vec<u8>>,
    pub q_tab: Vec<u8>,
    pub p_tab: Option<Vec<u8>>,
    pub d_tab: Option<Vec<u8>>,
}

pub fn fqz_decode_single_param(src: &mut &[u8]) -> io::Result<Parameter> {
    let context = read_u16_le(src)?;
    let flags = read_u8(src).map(Flags::from)?;
    let max_symbol = read_max_symbol(src)?;

    let (q_bits, q_shift) = read_u4x2(src)?;
    let (q_loc, s_loc) = read_u4x2(src)?;
    let (p_loc, d_loc) = read_u4x2(src)?;

    let q_map = if flags.contains(Flags::HAVE_QMAP) {
        read_quality_map(src, max_symbol).map(Some)?
    } else {
        None
    };

    let q_tab = if flags.contains(Flags::HAVE_QTAB) {
        read_array(src, 256)?
    } else {
        build_default_qualities_table()
    };

    let p_tab = if flags.contains(Flags::HAVE_PTAB) {
        read_array(src, 1024).map(Some)?
    } else {
        None
    };

    let d_tab = if flags.contains(Flags::HAVE_DTAB) {
        read_array(src, 256).map(Some)?
    } else {
        None
    };

    Ok(Parameter {
        context,
        flags,
        symbol_count: max_to_count(max_symbol),
        q_bits,
        q_shift,
        q_loc,
        s_loc,
        p_loc,
        d_loc,
        q_map,
        q_tab,
        p_tab,
        d_tab,
    })
}

fn read_u4x2(src: &mut &[u8]) -> io::Result<(u8, u8)> {
    let n = read_u8(src)?;
    Ok((n >> 4, n & 0x0f))
}

fn read_max_symbol(src: &mut &[u8]) -> io::Result<u8> {
    read_u8(src)
}

fn read_quality_map(src: &mut &[u8], max_symbol: u8) -> io::Result<Vec<u8>> {
    let len = usize::from(max_symbol);

    let (buf, rest) = src
        .split_at_checked(len)
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Ok(buf.to_vec())
}

fn build_default_qualities_table() -> Vec<u8> {
    (0..=u8::MAX).collect()
}

fn max_to_count(n: u8) -> NonZero<usize> {
    let m = usize::from(n) + 1;
    // SAFETY: `m > 0`.
    NonZero::new(m).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_u4x2() -> io::Result<()> {
        assert_eq!(read_u4x2(&mut &[0x0f][..])?, (0x00, 0x0f));
        Ok(())
    }
}
