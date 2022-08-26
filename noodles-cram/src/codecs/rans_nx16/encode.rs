mod order_0;
mod order_1;

use std::io::{self, Write};

use byteorder::WriteBytesExt;

use super::Flags;
use crate::writer::num::write_uint7;

#[allow(dead_code)]
pub fn rans_encode_nx16(flags: Flags, src: &[u8]) -> io::Result<Vec<u8>> {
    let mut dst = Vec::new();

    dst.write_u8(u8::from(flags))?;

    if !flags.contains(Flags::NO_SIZE) {
        let n =
            u32::try_from(src.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_uint7(&mut dst, n)?;
    }

    let n = if flags.contains(Flags::N32) { 32 } else { 4 };

    if flags.contains(Flags::STRIPE) {
        todo!("rans_encode_stripe");
    }

    if flags.contains(Flags::PACK) {
        todo!("encode_pack_meta");
    }

    if flags.contains(Flags::RLE) {
        todo!("encode_rle_meta");
    }

    if flags.contains(Flags::CAT) {
        dst.write_all(src)?;
    } else if flags.contains(Flags::ORDER) {
        let (normalized_contexts, compressed_data) = order_1::encode(src, n)?;
        // bits = 12, no compression (0)
        dst.write_u8(12 << 4)?;
        order_1::write_contexts(&mut dst, &normalized_contexts)?;
        dst.write_all(&compressed_data)?;
    } else {
        let (normalized_frequencies, compressed_data) = order_0::encode(src, n)?;
        order_0::write_frequencies(&mut dst, &normalized_frequencies)?;
        dst.write_all(&compressed_data)?;
    }

    Ok(dst)
}

fn build_frequencies(src: &[u8]) -> Vec<u32> {
    let mut frequencies = vec![0; 256];

    for &b in src {
        let i = usize::from(b);
        frequencies[i] += 1;
    }

    frequencies
}

fn normalize_frequencies(frequencies: &[u32]) -> Vec<u32> {
    const SCALE: u32 = 4096;

    let mut sum = 0;
    let mut max = 0;
    let mut max_index = 0;

    for (i, &f) in frequencies.iter().enumerate() {
        if f >= max {
            max = f;
            max_index = i;
        }

        sum += f;
    }

    if sum == 0 {
        return vec![0; frequencies.len()];
    }

    let mut normalize_frequencies = vec![0; frequencies.len()];
    let mut normalized_sum = 0;

    for (&f, g) in frequencies.iter().zip(normalize_frequencies.iter_mut()) {
        let normalized_frequency = f * SCALE / sum;
        *g = normalized_frequency;
        normalized_sum += normalized_frequency;
    }

    if normalized_sum < SCALE {
        normalize_frequencies[max_index] += SCALE - normalized_sum;
    }

    normalize_frequencies
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
            writer.write_u8(sym as u8)?;

            if sym > 0 && frequencies[sym - 1] > 0 {
                rle = frequencies[sym + 1..]
                    .iter()
                    .position(|&f| f == 0)
                    .unwrap_or(255);

                writer.write_u8(rle as u8)?;
            }
        }
    }

    writer.write_u8(0x00)?;

    Ok(())
}

pub fn normalize<W>(writer: &mut W, mut r: u32, freq_i: u32, bits: u32) -> io::Result<u32>
where
    W: Write,
{
    while r >= ((1 << (31 - bits)) * freq_i) {
        writer.write_u8(((r >> 8) & 0xff) as u8)?;
        writer.write_u8((r & 0xff) as u8)?;
        r >>= 16;
    }

    Ok(r)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rans_encode_nx16_order_0() -> io::Result<()> {
        let actual = rans_encode_nx16(Flags::empty(), b"noodles")?;

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
    fn test_rans_encode_nx16_uncompressed() -> io::Result<()> {
        let actual = rans_encode_nx16(Flags::CAT, b"noodles")?;

        let expected = [
            0x20, // flags = CAT
            0x07, // uncompressed len = 7
            0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
