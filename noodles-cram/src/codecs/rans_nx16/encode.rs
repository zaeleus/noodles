mod order_0;
mod order_1;
pub mod pack;
mod rle;

use std::io::{self, Write};

use byteorder::WriteBytesExt;

use super::Flags;
use crate::io::writer::num::write_uint7;

pub fn encode(mut flags: Flags, src: &[u8]) -> io::Result<Vec<u8>> {
    let mut src = src.to_vec();
    let mut dst = Vec::new();

    dst.write_u8(u8::from(flags))?;

    if !flags.contains(Flags::NO_SIZE) {
        let n =
            u32::try_from(src.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_uint7(&mut dst, n)?;
    }

    let n = if flags.contains(Flags::N32) { 32 } else { 4 };

    if flags.contains(Flags::STRIPE) {
        let buf = rans_encode_stripe(&src, n)?;
        dst.extend(&buf);
        return Ok(dst);
    }

    let mut pack_header = None;

    if flags.contains(Flags::PACK) {
        match pack::encode(&src) {
            Ok((header, buf)) => {
                pack_header = Some(header);
                src = buf;
            }
            Err(e) if e.kind() == io::ErrorKind::InvalidInput => {
                flags.remove(Flags::PACK);
                dst[0] = u8::from(flags);
            }
            Err(e) => return Err(e),
        }
    }

    let mut rle_header = None;

    if flags.contains(Flags::RLE) {
        let (header, buf) = rle::encode(&src)?;
        rle_header = Some(header);
        src = buf;
    }

    if let Some(header) = pack_header {
        dst.write_all(&header)?;
    }

    if let Some(header) = rle_header {
        dst.write_all(&header)?;
    }

    if src.len() < n {
        flags.remove(Flags::ORDER);

        if src.is_empty() {
            flags.insert(Flags::CAT);
        }

        dst[0] = u8::from(flags);
    }

    if flags.contains(Flags::CAT) {
        dst.write_all(&src)?;
    } else if flags.contains(Flags::ORDER) {
        let (normalized_contexts, compressed_data) = order_1::encode(&src, n)?;
        // bits = 12, no compression (0)
        dst.write_u8(12 << 4)?;
        order_1::write_contexts(&mut dst, &normalized_contexts)?;
        dst.write_all(&compressed_data)?;
    } else {
        let (normalized_frequencies, compressed_data) = order_0::encode(&src, n)?;
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
    use std::cmp::Ordering;

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

    let mut normalized_frequencies = vec![0; frequencies.len()];
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
                    .unwrap_or(0);

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

fn rans_encode_stripe(src: &[u8], n: usize) -> io::Result<Vec<u8>> {
    let mut ulens = Vec::with_capacity(n);
    let mut t = Vec::with_capacity(n);

    for j in 0..n {
        let mut ulen = src.len() / n;

        if src.len() % n > j {
            ulen += 1;
        }

        let chunk = vec![0; ulen];

        ulens.push(ulen);
        t.push(chunk);
    }

    let mut x = 0;
    let mut i = 0;

    while i < src.len() {
        for j in 0..n {
            if x < ulens[j] {
                t[j][x] = src[i + j];
            }
        }

        x += 1;
        i += n;
    }

    let mut chunks = vec![Vec::new(); n];

    for (chunk, s) in chunks.iter_mut().zip(t.iter()) {
        *chunk = encode(Flags::empty(), s)?;
    }

    let mut dst = Vec::new();

    dst.write_u8(n as u8)?;

    for chunk in &chunks {
        let clen = chunk.len() as u32;
        write_uint7(&mut dst, clen)?;
    }

    for chunk in &chunks {
        dst.write_all(chunk)?;
    }

    Ok(dst)
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
            0x08, 0x07, 0x04, 0x19, 0x19, 0x19, 0x16, 0x00, 0x02, 0x6c, 0x6e, 0x00, 0x90, 0x00,
            0x90, 0x00, 0x00, 0x08, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x80, 0x00, 0x00, 0x00, 0x02, 0x65, 0x6f, 0x00, 0x90, 0x00, 0x90, 0x00, 0x00,
            0x08, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00,
            0x00, 0x00, 0x02, 0x6f, 0x73, 0x00, 0x90, 0x00, 0x90, 0x00, 0x00, 0x00, 0x01, 0x00,
            0x00, 0x08, 0x01, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x01,
            0x64, 0x00, 0xa0, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
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
