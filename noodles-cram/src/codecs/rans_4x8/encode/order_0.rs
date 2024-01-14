use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use crate::io::writer::num::write_itf8;

use super::{
    build_cumulative_frequencies, normalize, normalize_frequencies, update, BASE, LOWER_BOUND,
};

pub fn encode(src: &[u8]) -> io::Result<Vec<u8>> {
    use super::{write_header, Order};

    let frequencies = build_frequencies(src, BASE);

    let freq = normalize_frequencies(&frequencies);
    let cfreq = build_cumulative_frequencies(&freq);

    let mut buf = Vec::new();
    let mut states = [LOWER_BOUND; 4];

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

    for &state in &states {
        dst.write_u32::<LittleEndian>(state)?;
    }

    dst.extend(buf.iter().rev());

    let compressed_len = dst[9..].len();
    let mut writer = &mut dst[..9];
    write_header(&mut writer, Order::Zero, compressed_len, src.len())?;

    Ok(dst)
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
                    .unwrap_or(0);

                writer.write_u8(rle as u8)?;
            }
        }

        write_itf8(writer, f as i32)?;
    }

    writer.write_u8(0x00)?;

    Ok(())
}

fn build_frequencies(src: &[u8], bin_count: usize) -> Vec<u32> {
    let mut frequencies = vec![0; bin_count];

    for &b in src {
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
    fn test_build_frequencies() {
        let data = [0, 1, 1, 2, 2, 2];
        assert_eq!(build_frequencies(&data, 4), [1, 2, 3, 0]);
    }
}
