mod header;
mod order_0;
mod order_1;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use self::header::write_header;
use super::{Order, ALPHABET_SIZE, LOWER_BOUND, STATE_COUNT};

pub fn encode(order: Order, src: &[u8]) -> io::Result<Vec<u8>> {
    match order {
        Order::Zero => order_0::encode(src),
        Order::One => order_1::encode(src),
    }
}

fn write_states<W>(writer: &mut W, states: &[u32; STATE_COUNT]) -> io::Result<()>
where
    W: Write,
{
    for state in states {
        writer.write_u32::<LittleEndian>(*state)?;
    }

    Ok(())
}

fn normalize_frequencies(raw_frequencies: &[u32]) -> [u16; ALPHABET_SIZE] {
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

fn build_cumulative_frequencies(frequencies: &[u16]) -> Vec<u16> {
    let mut cumulative_frequencies = vec![0; frequencies.len()];

    for i in 0..frequencies.len() - 1 {
        cumulative_frequencies[i + 1] = cumulative_frequencies[i] + frequencies[i];
    }

    cumulative_frequencies
}

fn normalize<W>(writer: &mut W, mut x: u32, freq_i: u16) -> io::Result<u32>
where
    W: Write,
{
    while x >= (LOWER_BOUND >> 4) * u32::from(freq_i) {
        let b = (x & 0xff) as u8;
        writer.write_u8(b)?;
        x >>= 8;
    }

    Ok(x)
}

fn update(x: u32, freq_i: u16, cfreq_i: u16) -> u32 {
    let (q, r) = (x / u32::from(freq_i), x % u32::from(freq_i));
    (q << 12) + r + u32::from(cfreq_i)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codecs::rans_4x8::ALPHABET_SIZE;

    #[test]
    fn test_encode_with_order_0() -> io::Result<()> {
        let data = b"noodles";
        let actual = encode(Order::Zero, data)?;

        let expected = [
            0x00, 0x25, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x64, 0x82, 0x49, 0x65, 0x00,
            0x82, 0x49, 0x6c, 0x82, 0x49, 0x6e, 0x82, 0x49, 0x6f, 0x00, 0x84, 0x92, 0x73, 0x82,
            0x49, 0x00, 0xe2, 0x06, 0x83, 0x18, 0x74, 0x7b, 0x41, 0x0c, 0x2b, 0xa9, 0x41, 0x0c,
            0x25, 0x31, 0x80, 0x03,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_with_order_1() -> io::Result<()> {
        let data = b"noodles";
        let actual = encode(Order::One, data)?;

        let expected = [
            0x01, 0x3b, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x00, 0x64, 0x83, 0xff, 0x6e,
            0x83, 0xff, 0x6f, 0x00, 0x88, 0x01, 0x00, 0x64, 0x6c, 0x8f, 0xff, 0x00, 0x65, 0x00,
            0x73, 0x8f, 0xff, 0x00, 0x6c, 0x65, 0x8f, 0xff, 0x00, 0x6e, 0x6f, 0x8f, 0xff, 0x00,
            0x6f, 0x00, 0x64, 0x87, 0xff, 0x6f, 0x88, 0x00, 0x00, 0x00, 0x07, 0x84, 0x00, 0x02,
            0x00, 0xe8, 0xff, 0x00, 0x00, 0xe8, 0xff, 0x00, 0x10, 0xe0, 0x00, 0x02,
        ];

        assert_eq!(actual, expected);

        Ok(())
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
