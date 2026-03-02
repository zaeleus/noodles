mod flags;

use std::{io, num::NonZero};

pub use self::flags::Flags;
use self::flags::read_flags;
use super::read_array;
use crate::io::reader::num::{read_u8, read_u16_le};

// § 6.2 "FQZComp Data Stream" (2023-03-15).
const QUALITIES_TABLE_SIZE: usize = 256;
const POSITIONS_TABLE_SIZE: usize = 1024;
const DELTAS_TABLE_SIZE: usize = 256;

pub struct Parameter {
    pub context: u16,
    flags: Flags,
    symbol_count: NonZero<usize>,
    pub q_bits: u8,
    pub q_shift: u8,
    pub q_loc: u8,
    pub s_loc: u8,
    pub p_loc: u8,
    pub d_loc: u8,
    quality_map: Option<Vec<u8>>,
    qualities_table: Vec<u8>,
    positions_table: Option<Vec<u8>>,
    deltas_table: Option<Vec<u8>>,
}

impl Parameter {
    pub fn flags(&self) -> Flags {
        self.flags
    }

    pub fn symbol_count(&self) -> NonZero<usize> {
        self.symbol_count
    }

    pub fn quality_map(&self) -> Option<&[u8]> {
        self.quality_map.as_deref()
    }

    pub fn qualities_table(&self) -> &[u8] {
        &self.qualities_table
    }

    pub fn positions_table(&self) -> Option<&[u8]> {
        self.positions_table.as_deref()
    }

    pub fn deltas_table(&self) -> Option<&[u8]> {
        self.deltas_table.as_deref()
    }
}

pub fn fqz_decode_single_param(src: &mut &[u8]) -> io::Result<Parameter> {
    let context = read_u16_le(src)?;
    let flags = read_flags(src)?;
    let max_symbol = read_max_symbol(src)?;

    let (q_bits, q_shift) = read_u4x2(src)?;
    let (q_loc, s_loc) = read_u4x2(src)?;
    let (p_loc, d_loc) = read_u4x2(src)?;

    let quality_map = if flags.has_quality_map() {
        read_quality_map(src, max_symbol).map(Some)?
    } else {
        None
    };

    let qualities_table = if flags.has_qualities_table() {
        read_qualities_table(src)?
    } else {
        build_default_qualities_table()
    };

    let positions_table = if flags.has_positions_table() {
        read_positions_table(src).map(Some)?
    } else {
        None
    };

    let deltas_table = if flags.has_deltas_table() {
        read_deltas_table(src).map(Some)?
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
        quality_map,
        qualities_table,
        positions_table,
        deltas_table,
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
    let len = usize::from(max_symbol) + 1;

    let (buf, rest) = src
        .split_at_checked(len)
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Ok(buf.to_vec())
}

fn read_qualities_table(src: &mut &[u8]) -> io::Result<Vec<u8>> {
    read_array(src, QUALITIES_TABLE_SIZE)
}

fn build_default_qualities_table() -> Vec<u8> {
    (0..=u8::MAX).collect()
}

fn read_positions_table(src: &mut &[u8]) -> io::Result<Vec<u8>> {
    read_array(src, POSITIONS_TABLE_SIZE)
}

fn read_deltas_table(src: &mut &[u8]) -> io::Result<Vec<u8>> {
    read_array(src, DELTAS_TABLE_SIZE)
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
