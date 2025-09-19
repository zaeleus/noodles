mod order_0;
mod order_1;
pub mod pack;
mod rle;

use std::{
    io::{self, Read},
    mem,
    num::NonZero,
};

use super::Flags;
use crate::io::reader::num::{read_u8, read_uint7};

pub fn decode<R>(reader: &mut R, mut len: usize) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let flags = read_flags(reader)?;

    let state_count = flags.state_count();

    if flags.has_uncompressed_size() {
        len = read_uint7(reader).and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;
    }

    if flags.is_striped() {
        return rans_decode_stripe(reader, len);
    }

    let mut p = None;
    let mut n_sym = None;
    let pack_len = len;

    if flags.is_bit_packed() {
        let (q, n, new_len) = decode_pack_meta(reader)?;
        p = Some(q);
        n_sym = Some(n);
        len = new_len;
    }

    let mut l = None;
    let mut rle_meta = None;
    let rle_len = len;

    if flags.is_rle() {
        let (m, meta, new_len) = rle::decode_rle_meta(reader, state_count)?;
        l = Some(m);
        rle_meta = Some(meta);
        len = new_len;
    }

    let mut data = vec![0; len];

    if flags.is_uncompressed() {
        reader.read_exact(&mut data)?;
    } else if flags.order() == 0 {
        order_0::decode(reader, &mut data, state_count)?;
    } else {
        order_1::decode(reader, &mut data, state_count)?;
    };

    if flags.is_rle() {
        let l = l.unwrap();
        let mut rle_meta = rle_meta.unwrap();
        data = rle::decode(&data, &l, &mut rle_meta, rle_len)?;
    }

    if flags.is_bit_packed() {
        let p = p.unwrap();
        let n_sym = n_sym.unwrap();
        data = pack::decode(&data, &p, n_sym, pack_len)?;
    }

    Ok(data)
}

fn read_flags<R>(reader: &mut R) -> io::Result<Flags>
where
    R: Read,
{
    read_u8(reader).map(Flags::from)
}

fn read_alphabet<R>(reader: &mut R) -> io::Result<[bool; 256]>
where
    R: Read,
{
    let mut alphabet = [false; 256];

    let mut sym = read_u8(reader)?;
    let mut last_sym = sym;
    let mut rle = 0;

    loop {
        alphabet[usize::from(sym)] = true;

        if rle > 0 {
            rle -= 1;
            sym += 1;
        } else {
            sym = read_u8(reader)?;

            if last_sym < 255 && sym == last_sym + 1 {
                rle = read_u8(reader)?;
            }
        }

        last_sym = sym;

        if sym == 0 {
            break;
        }
    }

    Ok(alphabet)
}

fn rans_get_cumulative_freq_nx16(r: u32, bits: u32) -> u32 {
    r & ((1 << bits) - 1)
}

fn rans_get_symbol_from_freq(cumulative_freqs: &[u32], freq: u32) -> u8 {
    let mut sym = 0;

    while sym < 255 && freq >= cumulative_freqs[(sym + 1) as usize] {
        sym += 1;
    }

    sym
}

fn rans_advance_step_nx16(r: u32, c: u32, f: u32, bits: u32) -> u32 {
    f * (r >> bits) + (r & ((1 << bits) - 1)) - c
}

pub fn rans_renorm_nx16<R>(reader: &mut R, mut r: u32) -> io::Result<u32>
where
    R: Read,
{
    if r < (1 << 15) {
        r = (r << 16) + read_u16_le(reader).map(u32::from)?;
    }

    Ok(r)
}

fn rans_decode_stripe<R>(reader: &mut R, len: usize) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let x = read_u8(reader).map(usize::from)?;
    let mut clens = Vec::with_capacity(x);

    for _ in 0..x {
        let clen = read_uint7(reader).and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        clens.push(clen);
    }

    let mut ulens = Vec::with_capacity(x);
    let mut t = Vec::with_capacity(x);

    for j in 0..x {
        let mut ulen = len / x;

        if len % x > j {
            ulen += 1;
        }

        let chunk = decode(reader, ulen)?;

        ulens.push(ulen);
        t.push(chunk);
    }

    let mut dst = vec![0; len];

    for j in 0..x {
        for i in 0..ulens[j] {
            dst[i * x + j] = t[j][i];
        }
    }

    Ok(dst)
}

pub fn decode_pack_meta<R>(reader: &mut R) -> io::Result<(Vec<u8>, NonZero<usize>, usize)>
where
    R: Read,
{
    // n_sym
    let symbol_count = read_u8(reader).and_then(|n| {
        NonZero::try_from(usize::from(n)).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut p = vec![0; symbol_count.get()];
    reader.read_exact(&mut p)?;

    let len = read_uint7(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    Ok((p, symbol_count, len))
}

fn read_u16_le<R>(reader: &mut R) -> io::Result<u16>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u16>()];
    reader.read_exact(&mut buf)?;
    Ok(u16::from_le_bytes(buf))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_order_0() -> io::Result<()> {
        let data = [
            0x00, // flags = {empty}
            0x07, // uncompressed len = 7
            0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x01, 0x01, 0x01, 0x01, 0x03,
            0x01, 0x00, 0x26, 0x20, 0x00, 0x00, 0xb8, 0x0a, 0x00, 0x00, 0xd8, 0x0a, 0x00, 0x00,
            0x00, 0x04, 0x00,
        ];
        let mut reader = &data[..];

        assert_eq!(decode(&mut reader, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_decode_order_1() -> io::Result<()> {
        let data = [
            0x01, // flags = ORDER
            0x07, // uncompressed len = 7
            0xc0, 0x00, 0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x00, 0x00, 0x01,
            0x00, 0x01, 0x01, 0x02, 0x00, 0x00, 0x00, 0x02, 0x01, 0x00, 0x02, 0x00, 0x05, 0x01,
            0x00, 0x01, 0x01, 0x00, 0x03, 0x00, 0x04, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
            0x02, 0x01, 0x00, 0x00, 0x01, 0x00, 0x05, 0x00, 0x04, 0x02, 0x00, 0x00, 0x08, 0x01,
            0x00, 0x00, 0x08, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00,
        ];

        let mut reader = &data[..];
        assert_eq!(decode(&mut reader, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_decode_stripe() -> io::Result<()> {
        let data = [
            0x08, // flags = STRIPE
            0x07, // uncompressed len = 7
            0x04, 0x17, 0x17, 0x17, 0x15, 0x00, 0x02, 0x6c, 0x6e, 0x00, 0x01, 0x01, 0x00, 0x08,
            0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x02, 0x65, 0x6f, 0x00, 0x01, 0x01, 0x00, 0x08, 0x01, 0x00, 0x00, 0x00, 0x01,
            0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x02, 0x6f, 0x73, 0x00,
            0x01, 0x01, 0x00, 0x00, 0x01, 0x00, 0x00, 0x08, 0x01, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x80, 0x00, 0x00, 0x00, 0x01, 0x64, 0x00, 0x01, 0x00, 0x80, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x02, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x22, 0x00, 0x81, 0x11, 0x01, 0x7f, 0x00,
        ];

        let mut reader = &data[..];
        assert_eq!(decode(&mut reader, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_decode_uncompressed() -> io::Result<()> {
        let data = [
            0x20, // flags = CAT
            0x07, // uncompressed len = 7
            0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73,
        ];

        let mut reader = &data[..];
        assert_eq!(decode(&mut reader, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_decode_rle() -> io::Result<()> {
        let data = [
            0x40, // flags = RLE
            0x0d, // uncompressed len = 13
            0x06, 0x06, 0x17, 0x01, 0x07, 0x6f, 0x00, 0x02, 0x01, 0x01, 0x00, 0x00, 0x01, 0x00,
            0x00, 0x0c, 0x02, 0x00, 0x00, 0x08, 0x02, 0x00, 0x00, 0x80, 0x00, 0x00, 0x64, 0x65,
            0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x03, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00,
            0x3a, 0x20, 0x00, 0x00, 0x7c, 0x20, 0x00, 0x00, 0x52, 0x01, 0x00, 0x00, 0x08, 0x04,
            0x00,
        ];

        let mut reader = &data[..];
        assert_eq!(decode(&mut reader, 0)?, b"noooooooodles");

        Ok(())
    }

    #[test]
    fn test_decode_bit_packing_with_6_symbols() -> io::Result<()> {
        let data = [
            0x80, // flags = PACK
            0x07, // uncompressed len = 7
            0x06, 0x64, 0x65, 0x6c, 0x6e, 0x6f, 0x73, 0x04, 0x04, 0x05, 0x00, 0x12, 0x43, 0x00,
            0x01, 0x01, 0x01, 0x01, 0x00, 0x0c, 0x02, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x08,
            0x02, 0x00, 0x00, 0x04, 0x02, 0x00,
        ];

        let mut reader = &data[..];
        assert_eq!(decode(&mut reader, 0)?, b"noodles");

        Ok(())
    }
}
