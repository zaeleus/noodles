use std::{
    io::{self, Write},
    mem,
};

use byteorder::{LittleEndian, WriteBytesExt};

use crate::writer::num::write_itf8;

use super::{
    build_cumulative_frequencies, normalize, normalize_frequencies, update, BASE, LOWER_BOUND,
};

pub fn encode(data: &[u8]) -> io::Result<(Vec<u32>, Vec<u8>)> {
    let frequencies = build_frequencies(data, BASE);

    let freq = normalize_frequencies(&frequencies);
    let cfreq = build_cumulative_frequencies(&freq);

    let mut buf = Vec::new();
    let mut states = [LOWER_BOUND; 4];

    for (i, &sym) in data.iter().enumerate().rev() {
        let j = i % states.len();

        let mut x = states[j];
        let freq_i = freq[usize::from(sym)];
        let cfreq_i = cfreq[usize::from(sym)];

        x = normalize(&mut buf, x, freq_i)?;
        states[j] = update(x, freq_i, cfreq_i);
    }

    let mut writer = Vec::with_capacity(states.len() * mem::size_of::<u32>() + buf.len());

    for &state in &states {
        writer.write_u32::<LittleEndian>(state)?;
    }

    writer.extend(buf.iter().rev());

    Ok((freq, writer))
}

pub fn write_frequencies<W>(writer: &mut W, frequencies: &[u32]) -> io::Result<()>
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

fn build_frequencies(data: &[u8], bin_count: usize) -> Vec<u32> {
    let mut frequencies = vec![0; bin_count];

    for &b in data {
        let i = usize::from(b);
        frequencies[i] += 1;
    }

    frequencies
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode() -> io::Result<()> {
        // § 14.4.2 Frequency table, Order-0 encoding (2020-07-22)
        let data = b"abracadabra";
        let (actual_normalized_frequencies, actual_compressed_data) = encode(data)?;

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

    #[test]
    fn test_build_frequencies() {
        let data = [0, 1, 1, 2, 2, 2];
        assert_eq!(build_frequencies(&data, 4), [1, 2, 3, 0]);
    }
}
