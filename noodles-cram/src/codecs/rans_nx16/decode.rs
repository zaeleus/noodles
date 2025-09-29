pub mod bit_pack;
mod order_0;
mod order_1;
mod rle;
mod stripe;

use std::io::{self, Read};

use super::{ALPHABET_SIZE, Flags};
use crate::io::reader::num::{read_u8, read_u32_le, read_uint7_as};

pub fn decode(mut src: &[u8], mut len: usize) -> io::Result<Vec<u8>> {
    let flags = read_flags(&mut src)?;

    let state_count = flags.state_count();

    if flags.has_uncompressed_size() {
        len = read_uint7_as(&mut src)?;
    }

    if flags.is_striped() {
        return stripe::decode(&mut src, len);
    }

    let mut bit_pack_context = None;

    if flags.is_bit_packed() {
        let (ctx, new_len) = bit_pack::read_context(&mut src, len)?;
        bit_pack_context = Some(ctx);
        len = new_len;
    }

    let mut rle_context = None;

    if flags.is_rle() {
        let (ctx, new_len) = rle::read_context(&mut src, state_count, len)?;
        rle_context = Some(ctx);
        len = new_len;
    }

    let mut dst = vec![0; len];

    if flags.is_uncompressed() {
        src.read_exact(&mut dst)?;
    } else if flags.order() == 0 {
        order_0::decode(&mut src, &mut dst, state_count)?;
    } else {
        order_1::decode(&mut src, &mut dst, state_count)?;
    };

    if let Some(ctx) = rle_context {
        dst = rle::decode(&dst, &ctx)?;
    }

    if let Some(ctx) = bit_pack_context {
        dst = bit_pack::decode(&dst, &ctx)?;
    }

    Ok(dst)
}

fn read_flags(src: &mut &[u8]) -> io::Result<Flags> {
    read_u8(src).map(Flags::from)
}

fn read_alphabet(src: &mut &[u8]) -> io::Result<[bool; ALPHABET_SIZE]> {
    let mut alphabet = [false; ALPHABET_SIZE];

    let mut sym = read_u8(src)?;
    let mut last_sym = sym;
    let mut rle = 0;

    loop {
        alphabet[usize::from(sym)] = true;

        if rle > 0 {
            rle -= 1;
            sym += 1;
        } else {
            sym = read_u8(src)?;

            if last_sym < 255 && sym == last_sym + 1 {
                rle = read_u8(src)?;
            }
        }

        last_sym = sym;

        if sym == 0 {
            break;
        }
    }

    Ok(alphabet)
}

fn read_states(src: &mut &[u8], state_count: usize) -> io::Result<Vec<u32>> {
    (0..state_count).map(|_| read_u32_le(src)).collect()
}

fn state_cumulative_frequency(s: u32, bits: u32) -> u32 {
    s & ((1 << bits) - 1)
}

fn cumulative_frequencies_symbol(cumulative_frequencies: &[u32], frequency: u32) -> u8 {
    let mut sym = 0;

    while sym < 255 && frequency >= cumulative_frequencies[usize::from(sym + 1)] {
        sym += 1;
    }

    sym
}

fn state_step(s: u32, f: u32, g: u32, bits: u32) -> u32 {
    f * (s >> bits) + (s & ((1 << bits) - 1)) - g
}

fn state_renormalize(mut s: u32, src: &mut &[u8]) -> io::Result<u32> {
    if s < (1 << 15) {
        s = (s << 16) + read_u16_le(src).map(u32::from)?;
    }

    Ok(s)
}

fn read_u16_le(src: &mut &[u8]) -> io::Result<u16> {
    let (buf, rest) = src
        .split_first_chunk()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Ok(u16::from_le_bytes(*buf))
}

fn split_off<'a>(src: &mut &'a [u8], len: usize) -> io::Result<&'a [u8]> {
    let (buf, rest) = src
        .split_at_checked(len)
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_order_0() -> io::Result<()> {
        let src = [
            0x00, // flags = {empty}
            0x07, // uncompressed len = 7
            0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x01, 0x01, 0x01, 0x01, 0x03,
            0x01, 0x00, 0x26, 0x20, 0x00, 0x00, 0xb8, 0x0a, 0x00, 0x00, 0xd8, 0x0a, 0x00, 0x00,
            0x00, 0x04, 0x00,
        ];

        assert_eq!(decode(&src, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_decode_order_1() -> io::Result<()> {
        let src = [
            1, 77, 160, 0, 100, 101, 0, 108, 110, 111, 0, 115, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0,
            0, 15, 0, 0, 1, 0, 2, 0, 1, 15, 0, 2, 1, 0, 1, 1, 15, 0, 2, 0, 3, 15, 1, 0, 0, 0, 0, 1,
            0, 2, 15, 0, 0, 0, 5, 16, 128, 114, 96, 0, 128, 139, 95, 0, 192, 176, 96, 0, 64, 73,
            57, 0,
        ];

        assert_eq!(
            decode(&src, 0)?,
            b"nnnnnnnnnnnnooooooooooooooooddddddddddddddllllllllllllllleeeeeeeeeessssssssss"
        );

        Ok(())
    }

    #[test]
    fn test_decode_stripe() -> io::Result<()> {
        let src = [
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

        assert_eq!(decode(&src, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_decode_uncompressed() -> io::Result<()> {
        let src = [
            0x20, // flags = CAT
            0x07, // uncompressed len = 7
            0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73,
        ];

        assert_eq!(decode(&src, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_decode_rle() -> io::Result<()> {
        let src = [
            0x40, // flags = RLE
            0x0d, // uncompressed len = 13
            0x06, 0x06, 0x17, 0x01, 0x07, 0x6f, 0x00, 0x02, 0x01, 0x01, 0x00, 0x00, 0x01, 0x00,
            0x00, 0x0c, 0x02, 0x00, 0x00, 0x08, 0x02, 0x00, 0x00, 0x80, 0x00, 0x00, 0x64, 0x65,
            0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x03, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00,
            0x3a, 0x20, 0x00, 0x00, 0x7c, 0x20, 0x00, 0x00, 0x52, 0x01, 0x00, 0x00, 0x08, 0x04,
            0x00,
        ];

        assert_eq!(decode(&src, 0)?, b"noooooooodles");

        Ok(())
    }

    #[test]
    fn test_decode_bit_packing_with_6_symbols() -> io::Result<()> {
        let src = [
            0x80, // flags = PACK
            0x07, // uncompressed len = 7
            0x06, 0x64, 0x65, 0x6c, 0x6e, 0x6f, 0x73, 0x04, 0x04, 0x05, 0x00, 0x12, 0x43, 0x00,
            0x01, 0x01, 0x01, 0x01, 0x00, 0x0c, 0x02, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x08,
            0x02, 0x00, 0x00, 0x04, 0x02, 0x00,
        ];

        assert_eq!(decode(&src, 0)?, b"noodles");

        Ok(())
    }
}
