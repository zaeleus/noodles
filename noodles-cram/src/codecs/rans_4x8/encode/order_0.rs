use std::io::{self, Write};

use byteorder::WriteBytesExt;

use super::{normalize, update, write_header, write_states};
use crate::{
    codecs::rans_4x8::{Order, ALPHABET_SIZE, LOWER_BOUND, STATE_COUNT},
    io::writer::num::write_itf8,
};

pub fn encode(src: &[u8]) -> io::Result<Vec<u8>> {
    let raw_frequencies = build_raw_frequencies(src);

    let freq = normalize_frequencies(&raw_frequencies);
    let cfreq = build_cumulative_frequencies(&freq);

    let mut buf = Vec::new();
    let mut states = [LOWER_BOUND; STATE_COUNT];

    for (i, &sym) in src.iter().enumerate().rev() {
        let j = i % states.len();

        let mut x = states[j];
        let freq_i = freq[usize::from(sym)];
        let cfreq_i = cfreq[usize::from(sym)];

        x = normalize(&mut buf, x, freq_i)?;
        states[j] = update(x, freq_i, cfreq_i);
    }

    let mut dst = vec![0; 9];

    write_frequencies(&mut dst, &freq)?;
    write_states(&mut dst, &states)?;
    dst.extend(buf.iter().rev());

    let compressed_len = dst[9..].len();
    let mut writer = &mut dst[..9];
    write_header(&mut writer, Order::Zero, compressed_len, src.len())?;

    Ok(dst)
}

pub fn write_frequencies<W>(writer: &mut W, frequencies: &[u16]) -> io::Result<()>
where
    W: Write,
{
    assert_eq!(frequencies.len(), ALPHABET_SIZE);

    const NUL: u8 = 0x00;

    let mut iter = frequencies.iter().enumerate();
    let mut prev_sym = 0;

    while let Some((sym, &f)) = iter.next() {
        if f == 0 {
            continue;
        }

        // SAFETY: `sym <= ALPHABET_SIZE`.
        writer.write_u8(sym as u8)?;

        if sym > 0 && sym - 1 == prev_sym {
            let i = sym + 1;
            let len = frequencies[i..].iter().position(|&g| g == 0).unwrap_or(0);

            // SAFETY: `len < ALPHABET_SIZE`.
            writer.write_u8(len as u8)?;

            write_itf8(writer, i32::from(f))?;

            for (sym, &g) in iter.by_ref().take(len) {
                write_itf8(writer, i32::from(g))?;
                prev_sym = sym;
            }

            continue;
        }

        write_itf8(writer, i32::from(f))?;

        prev_sym = sym;
    }

    writer.write_u8(NUL)?;

    Ok(())
}

fn build_raw_frequencies(src: &[u8]) -> [u32; ALPHABET_SIZE] {
    let mut frequencies = [0; ALPHABET_SIZE];

    for &b in src {
        let i = usize::from(b);
        frequencies[i] += 1;
    }

    frequencies
}

pub(super) fn normalize_frequencies(raw_frequencies: &[u32]) -> [u16; ALPHABET_SIZE] {
    use std::cmp::Ordering;

    // § 2.1 "Frequency table" (2023-03-15): "The total sum of symbol frequencies are normalised to
    // add up to 4095."
    const SCALING_FACTOR: u16 = 4095;

    let (max_index, sum) = describe_frequencies(raw_frequencies);

    if sum == 0 {
        return [0; ALPHABET_SIZE];
    }

    let mut normalized_frequencies = [0; ALPHABET_SIZE];
    let mut normalized_sum = 0;

    for (i, &f) in raw_frequencies.iter().enumerate() {
        if f == 0 {
            continue;
        }

        let mut normalized_frequency = f * u32::from(SCALING_FACTOR) / sum;

        if normalized_frequency == 0 {
            normalized_frequency = 1;
        }

        // SAFETY: `normalized_frequency <= SCALING_FACTOR`.
        normalized_frequencies[i] = normalized_frequency as u16;
        normalized_sum += normalized_frequencies[i];
    }

    match normalized_sum.cmp(&SCALING_FACTOR) {
        Ordering::Less => normalized_frequencies[max_index] += SCALING_FACTOR - normalized_sum,
        Ordering::Equal => {}
        Ordering::Greater => normalized_frequencies[max_index] -= normalized_sum - SCALING_FACTOR,
    }

    normalized_frequencies
}

fn describe_frequencies(raw_frequencies: &[u32]) -> (usize, u32) {
    let mut max = u32::MIN;
    let mut max_index = 0;
    let mut sum = 0;

    for (i, &f) in raw_frequencies.iter().enumerate() {
        if f >= max {
            max = f;
            max_index = i;
        }

        sum += f;
    }

    (max_index, sum)
}

pub(super) fn build_cumulative_frequencies(frequencies: &[u16]) -> Vec<u16> {
    let mut cumulative_frequencies = vec![0; frequencies.len()];

    for i in 0..frequencies.len() - 1 {
        cumulative_frequencies[i + 1] = cumulative_frequencies[i] + frequencies[i];
    }

    cumulative_frequencies
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode() -> io::Result<()> {
        // § 14.4.2 Frequency table, Order-0 encoding (2020-07-22)
        let data = b"abracadabra";
        let actual = encode(data)?;

        let expected = [
            // header
            0x00, // order = 0
            0x1f, 0x00, 0x00, 0x00, // compressed_len = 31
            0x0b, 0x00, 0x00, 0x00, // uncompressed_len = 11
            // frequency table
            0x61, // sym = 'a'
            0x87, 0x47, // f['a'] = 1863
            0x62, // sym = 'b'
            0x02, // rle = 2
            0x82, 0xe8, // f['b'] = 744
            0x81, 0x74, // f['c'] = 372
            0x81, 0x74, // f['d'] = 372
            0x72, // 'r'
            0x82, 0xe8, // f['r'] = 744
            0x00, // end
            // compressed blob
            0xd2, 0x02, 0xa4, 0x42, // states[0]
            0x0d, 0x3a, 0x52, 0x21, // states[1]
            0xd0, 0xfe, 0xa1, 0x42, // states[2]
            0x40, 0xa6, 0x6a, 0x02, // states[3]
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_build_raw_frequencies() {
        // § 2.1.1 "Frequency table: Order-0 encoding" (2023-03-15)
        let src = b"abracadabra";

        let mut expected = [0; ALPHABET_SIZE];
        expected[usize::from(b'a')] = 5;
        expected[usize::from(b'b')] = 2;
        expected[usize::from(b'c')] = 1;
        expected[usize::from(b'd')] = 1;
        expected[usize::from(b'r')] = 2;

        assert_eq!(build_raw_frequencies(src), expected);
    }

    #[test]
    fn test_normalize_frequencies() {
        let mut raw_frequencies = [0; ALPHABET_SIZE];
        raw_frequencies[usize::from(b'a')] = 5;
        raw_frequencies[usize::from(b'b')] = 2;
        raw_frequencies[usize::from(b'c')] = 1;
        raw_frequencies[usize::from(b'd')] = 1;
        raw_frequencies[usize::from(b'r')] = 2;

        let actual = normalize_frequencies(&raw_frequencies);

        let mut expected = [0; ALPHABET_SIZE];
        expected[usize::from(b'a')] = 1863;
        expected[usize::from(b'b')] = 744;
        expected[usize::from(b'c')] = 372;
        expected[usize::from(b'd')] = 372;
        expected[usize::from(b'r')] = 744;

        assert_eq!(actual, expected);

        let raw_frequencies = [0; ALPHABET_SIZE];
        assert_eq!(normalize_frequencies(&raw_frequencies), [0; ALPHABET_SIZE]);
    }

    #[test]
    fn test_build_cumulative_frequencies() {
        let frequencies = [682, 1365, 2048, 0];

        assert_eq!(
            build_cumulative_frequencies(&frequencies),
            [0, 682, 2047, 4095]
        );
    }
}
