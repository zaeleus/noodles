pub mod bit_pack;
mod order_0;
mod order_1;
mod rle;
mod stripe;

use std::io::{self, Write};

use super::{ALPHABET_SIZE, Flags};
use crate::io::writer::num::{write_u8, write_u32_le, write_uint7};

// § 3 "rANS Nx16" (2023-03-15): "The lower-bound and initial encoder state _L_ is [...] 0x8000."
const LOWER_BOUND: u32 = 0x8000;

pub fn encode(mut flags: Flags, src: &[u8]) -> io::Result<Vec<u8>> {
    let mut src = src.to_vec();
    let mut dst = Vec::new();

    write_flags(&mut dst, flags)?;

    if flags.has_uncompressed_size() {
        let n =
            u32::try_from(src.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_uint7(&mut dst, n)?;
    }

    let state_count = flags.state_count();

    if flags.is_striped() {
        let buf = stripe::encode(&src)?;
        dst.extend(&buf);
        return Ok(dst);
    }

    if flags.is_bit_packed() {
        match bit_pack::build_context(&src) {
            Ok(ctx) => {
                src = bit_pack::encode(&src, &ctx);
                bit_pack::write_context(&mut dst, &ctx, src.len())?;
            }
            Err(
                bit_pack::context::BuildContextError::EmptyAlphabet
                | bit_pack::context::BuildContextError::TooManySymbols(_),
            ) => {
                flags.remove(Flags::PACK);
                dst[0] = u8::from(flags);
            }
        }
    }

    if flags.is_rle() {
        let mut ctx =
            rle::build_context(&src).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        src = rle::encode(&src, &mut ctx)?;
        rle::write_context(&mut dst, &ctx, src.len())?;
    }

    if src.len() < state_count {
        flags.remove(Flags::ORDER);
        flags.insert(Flags::CAT);
        dst[0] = u8::from(flags);
    }

    if flags.is_uncompressed() {
        dst.write_all(&src)?;
    } else if flags.order() == 0 {
        let ctx = order_0::build_context(&src);
        let compressed_data = order_0::encode(&src, &ctx, state_count)?;
        order_0::write_context(&mut dst, &ctx)?;
        dst.write_all(&compressed_data)?;
    } else {
        let (normalized_contexts, compressed_data) = order_1::encode(&src, state_count)?;
        // bits = 12, no compression (0)
        write_u8(&mut dst, 12 << 4)?;
        order_1::write_contexts(&mut dst, &normalized_contexts)?;
        dst.write_all(&compressed_data)?;
    }

    Ok(dst)
}

fn write_flags(dst: &mut Vec<u8>, flags: Flags) -> io::Result<()> {
    write_u8(dst, u8::from(flags))
}

fn write_states(dst: &mut Vec<u8>, states: &[u32]) -> io::Result<()> {
    for &state in states {
        write_u32_le(dst, state)?;
    }

    Ok(())
}

fn build_frequencies(src: &[u8]) -> [u32; ALPHABET_SIZE] {
    let mut frequencies = [0; ALPHABET_SIZE];

    for &b in src {
        let i = usize::from(b);
        frequencies[i] += 1;
    }

    frequencies
}

fn normalize_frequencies(frequencies: &[u32]) -> [u32; ALPHABET_SIZE] {
    use std::cmp::Ordering;

    const SCALE: u32 = 4096;

    let (max_index, sum) = describe_frequencies(frequencies);

    if sum == 0 {
        return [0; ALPHABET_SIZE];
    }

    let mut normalized_frequencies = [0; ALPHABET_SIZE];
    let mut normalized_sum = 0;

    for (&f, g) in frequencies.iter().zip(normalized_frequencies.iter_mut()) {
        if f == 0 {
            continue;
        }

        let mut normalized_frequency = f * SCALE / sum;

        if normalized_frequency == 0 {
            normalized_frequency = 1;
        }

        *g = normalized_frequency;
        normalized_sum += normalized_frequency;
    }

    match normalized_sum.cmp(&SCALE) {
        Ordering::Less => normalized_frequencies[max_index] += SCALE - normalized_sum,
        Ordering::Equal => {}
        Ordering::Greater => normalized_frequencies[max_index] -= normalized_sum - SCALE,
    }

    normalized_frequencies
}

fn describe_frequencies(frequencies: &[u32]) -> (usize, u32) {
    let mut max = u32::MIN;
    let mut max_index = 0;
    let mut sum = 0;

    for (i, &f) in frequencies.iter().enumerate() {
        if f >= max {
            max = f;
            max_index = i;
        }

        sum += f;
    }

    (max_index, sum)
}

fn build_cumulative_frequencies(frequencies: &[u32]) -> Vec<u32> {
    let mut cumulative_frequencies = vec![0; frequencies.len()];

    for i in 0..frequencies.len() - 1 {
        cumulative_frequencies[i + 1] = cumulative_frequencies[i] + frequencies[i];
    }

    cumulative_frequencies
}

fn update(r: u32, c: u32, f: u32, bits: u32) -> u32 {
    ((r / f) << bits) + c + (r % f)
}

fn write_alphabet<W>(writer: &mut W, frequencies: &[u32]) -> io::Result<()>
where
    W: Write,
{
    let mut rle = 0;

    for (sym, &f) in frequencies.iter().enumerate() {
        if f == 0 {
            continue;
        }

        if rle > 0 {
            rle -= 1;
        } else {
            write_u8(writer, sym as u8)?;

            if sym > 0 && frequencies[sym - 1] > 0 {
                rle = frequencies[sym + 1..]
                    .iter()
                    .position(|&f| f == 0)
                    .unwrap_or(0);

                write_u8(writer, rle as u8)?;
            }
        }
    }

    write_u8(writer, 0x00)?;

    Ok(())
}

pub fn normalize<W>(writer: &mut W, mut r: u32, freq_i: u32, bits: u32) -> io::Result<u32>
where
    W: Write,
{
    while r >= ((1 << (31 - bits)) * freq_i) {
        write_u8(writer, ((r >> 8) & 0xff) as u8)?;
        write_u8(writer, (r & 0xff) as u8)?;
        r >>= 16;
    }

    Ok(r)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_order_0() -> io::Result<()> {
        let actual = encode(Flags::empty(), b"noodles")?;

        let expected = [
            0x00, // flags = {empty}
            0x07, // uncompressed len = 7
            0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x84, 0x49, 0x84, 0x49, 0x84,
            0x49, 0x84, 0x49, 0x89, 0x13, 0x84, 0x49, 0x1b, 0xa7, 0x18, 0x00, 0xe9, 0x4a, 0x0c,
            0x00, 0x31, 0x6d, 0x0c, 0x00, 0x08, 0x80, 0x03, 0x00,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_order_1() -> io::Result<()> {
        let actual = encode(Flags::ORDER, b"noodles")?;

        let expected = [
            0x01, 0x07, 0xc0, 0x00, 0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x00,
            0x00, 0x88, 0x00, 0x00, 0x01, 0x88, 0x00, 0x90, 0x00, 0x00, 0x00, 0x00, 0x02, 0xa0,
            0x00, 0x00, 0x02, 0x00, 0x05, 0xa0, 0x00, 0x00, 0x01, 0xa0, 0x00, 0x00, 0x03, 0x00,
            0x04, 0xa0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x90, 0x00, 0x00, 0x02, 0x90, 0x00, 0x00,
            0x00, 0xa0, 0x00, 0x00, 0x05, 0x00, 0x04, 0x02, 0x00, 0x00, 0x08, 0x01, 0x00, 0x00,
            0x08, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_stripe() -> io::Result<()> {
        let actual = encode(Flags::STRIPE, b"noodles")?;

        let expected = [
            0x08, // flags = STRIPE
            0x07, // uncompressed len = 7
            0x04, // chunk count = 4
            0x03, 0x03, 0x03, 0x02, // compressed sizes = [3, 3, 3, 2]
            0x30, 0x6e, 0x6c, // chunks[0]
            0x30, 0x6f, 0x65, // chunks[1]
            0x30, 0x6f, 0x73, // chunks[2]
            0x30, 0x64, // chunks[3]
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_uncompressed() -> io::Result<()> {
        let actual = encode(Flags::CAT, b"noodles")?;

        let expected = [
            0x20, // flags = CAT
            0x07, // uncompressed len = 7
            0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_rle() -> io::Result<()> {
        let actual = encode(Flags::CAT | Flags::RLE, b"noooooooodles")?;

        let expected = [
            0x60, // flags = CAT | RLE
            0x0d, // uncompressed length = 13
            0x07, 0x06, 0x01, 0x6f, 0x07, 0x6e, 0x6f, 0x64, 0x6c, 0x65, 0x73,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_pack() -> io::Result<()> {
        let actual = encode(Flags::PACK, b"noodles")?;

        let expected = [
            0x80, // flags = PACK
            0x07, // uncompressed len = 7
            0x06, 0x64, 0x65, 0x6c, 0x6e, 0x6f, 0x73, 0x04, 0x04, 0x05, 0x00, 0x12, 0x43, 0x00,
            0x88, 0x00, 0x88, 0x00, 0x88, 0x00, 0x88, 0x00, 0x00, 0x0c, 0x02, 0x00, 0x00, 0x00,
            0x02, 0x00, 0x00, 0x08, 0x02, 0x00, 0x00, 0x04, 0x02, 0x00,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
