use std::{
    io::{self, Write},
    mem,
};

use byteorder::{LittleEndian, WriteBytesExt};

use crate::num::write_itf8;

// Base `b`.
const BASE: usize = 256;

// Lower bound `L`.
const LOWER_BOUND: u32 = 0x800000;

fn build_frequencies(data: &[u8], bin_count: usize) -> Vec<u32> {
    let mut frequencies = vec![0; bin_count];

    for &b in data {
        let i = usize::from(b);
        frequencies[i] += 1;
    }

    frequencies
}

fn normalize_frequencies(frequencies: &[u32]) -> Vec<u32> {
    const SCALE: u32 = 4095;

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

    let mut normalized_sum = 0;
    let mut normalized_frequencies = vec![0; frequencies.len()];

    for (i, &f) in frequencies.iter().enumerate() {
        let normalized_frequency = f * SCALE / sum;
        normalized_frequencies[i] = normalized_frequency;
        normalized_sum += normalized_frequency;
    }

    // Because the calculation of `normalized_frequency` uses integer division (truncation), it's
    // possible that the sum of all the normalized frequencies is smaller than the scale value. In
    // this case, the difference is added to the last max value.
    if normalized_sum < SCALE {
        normalized_frequencies[max_index] += SCALE - normalized_sum;
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

#[allow(dead_code)]
fn rans_encode_0(data: &[u8]) -> io::Result<(Vec<u32>, Vec<u8>)> {
    let frequencies = build_frequencies(data, BASE);

    let freq = normalize_frequencies(&frequencies);
    let cfreq = build_cumulative_frequencies(&freq);

    let mut buf = Vec::new();
    let mut states = [LOWER_BOUND; 4];

    for (i, &sym) in data.iter().enumerate().rev() {
        let j = i % states.len();

        let mut x = states[j];
        let freq_i = freq[sym as usize];
        let cfreq_i = cfreq[sym as usize];

        // normalize
        while x >= (LOWER_BOUND >> 4) * freq_i {
            let b = (x & 0xff) as u8;
            buf.write_u8(b)?;
            x >>= 8;
        }

        // update
        states[j] = (x / freq_i) * 0x1000 + cfreq_i + (x % freq_i);
    }

    let mut writer = Vec::with_capacity(states.len() * mem::size_of::<u32>() + buf.len());

    for &state in &states {
        writer.write_u32::<LittleEndian>(state)?;
    }

    writer.extend(buf.iter().rev());

    Ok((freq, writer))
}

#[allow(dead_code)]
fn write_frequencies<W>(writer: &mut W, frequencies: &[u32]) -> io::Result<()>
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

        write_itf8(writer, f as i32)?;
    }

    writer.write_u8(0x00)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_frequencies() {
        let data = [0, 1, 1, 2, 2, 2];
        assert_eq!(build_frequencies(&data, 4), [1, 2, 3, 0]);
    }

    #[test]
    fn test_normalize_frequencies() {
        let frequencies = [1, 2, 3, 0];
        assert_eq!(normalize_frequencies(&frequencies), [682, 1365, 2048, 0]);
    }

    #[test]
    fn test_build_cumulative_frequencies() {
        let frequencies = [682, 1365, 2048, 0];

        assert_eq!(
            build_cumulative_frequencies(&frequencies),
            [0, 682, 2047, 4095]
        );
    }

    #[test]
    fn test_rans_encode_0() -> io::Result<()> {
        // § 14.4.2 Frequency table, Order-0 encoding (2020-07-22)
        let data = b"abracadabra";
        let (actual_normalized_frequencies, actual_compressed_data) = rans_encode_0(data)?;

        let mut expected_normalized_frequencies = [0; BASE];
        expected_normalized_frequencies[usize::from(b'a')] = 1863;
        expected_normalized_frequencies[usize::from(b'b')] = 744;
        expected_normalized_frequencies[usize::from(b'c')] = 372;
        expected_normalized_frequencies[usize::from(b'd')] = 372;
        expected_normalized_frequencies[usize::from(b'r')] = 744;

        assert_eq!(
            actual_normalized_frequencies,
            expected_normalized_frequencies
        );

        let expected_compressed_data = [
            0xd2, 0x02, 0xa4, 0x42, // states[0]
            0x0d, 0x3a, 0x52, 0x21, // states[1]
            0xd0, 0xfe, 0xa1, 0x42, // states[2]
            0x40, 0xa6, 0x6a, 0x02, // states[3]
        ];

        assert_eq!(actual_compressed_data, expected_compressed_data);

        Ok(())
    }

    #[test]
    fn test_write_frequencies() -> io::Result<()> {
        let data = b"abracadabra";
        let frequencies = build_frequencies(data, BASE);
        let normalized_frequencies = normalize_frequencies(&frequencies);

        let mut writer = Vec::new();
        write_frequencies(&mut writer, &normalized_frequencies)?;

        let expected = [
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
        ];

        assert_eq!(writer, expected);

        Ok(())
    }
}
