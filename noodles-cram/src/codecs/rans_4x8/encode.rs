mod order_0;
mod order_1;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use super::{Order, LOWER_BOUND};

pub fn encode(order: Order, src: &[u8]) -> io::Result<Vec<u8>> {
    match order {
        Order::Zero => order_0::encode(src),
        Order::One => order_1::encode(src),
    }
}

fn write_header<W>(
    writer: &mut W,
    order: Order,
    compressed_len: usize,
    uncompressed_len: usize,
) -> io::Result<()>
where
    W: Write,
{
    writer.write_u8(u8::from(order))?;

    let compressed_size = u32::try_from(compressed_len)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(compressed_size)?;

    let data_size = u32::try_from(uncompressed_len)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(data_size)?;

    Ok(())
}

fn normalize_frequencies(frequencies: &[u32]) -> Vec<u32> {
    use std::cmp::Ordering;

    // § 2.1 "Frequency table" (2023-03-15): "The total sum of symbol frequencies are normalised to
    // add up to 4095."
    const SCALING_FACTOR: u32 = 4095;

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

    let mut normalized_sum = 0;
    let mut normalized_frequencies = vec![0; frequencies.len()];

    for (i, &f) in frequencies.iter().enumerate() {
        if f == 0 {
            continue;
        }

        let mut normalized_frequency = f * SCALING_FACTOR / sum;

        if normalized_frequency == 0 {
            normalized_frequency = 1;
        }

        normalized_frequencies[i] = normalized_frequency;
        normalized_sum += normalized_frequency;
    }

    match normalized_sum.cmp(&SCALING_FACTOR) {
        Ordering::Less => normalized_frequencies[max_index] += SCALING_FACTOR - normalized_sum,
        Ordering::Equal => {}
        Ordering::Greater => normalized_frequencies[max_index] -= normalized_sum - SCALING_FACTOR,
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

fn normalize<W>(writer: &mut W, mut x: u32, freq_i: u32) -> io::Result<u32>
where
    W: Write,
{
    while x >= (LOWER_BOUND >> 4) * freq_i {
        let b = (x & 0xff) as u8;
        writer.write_u8(b)?;
        x >>= 8;
    }

    Ok(x)
}

fn update(x: u32, freq_i: u32, cfreq_i: u32) -> u32 {
    (x / freq_i) * 0x1000 + cfreq_i + (x % freq_i)
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_write_header() -> io::Result<()> {
        let mut writer = Vec::new();
        write_header(&mut writer, Order::One, 14930352, 9227465)?;

        let expected = [
            0x01, // order
            0xb0, 0xd1, 0xe3, 0x00, // compressed length
            0xc9, 0xcc, 0x8c, 0x00, // data length
        ];

        assert_eq!(writer, expected);

        Ok(())
    }

    #[test]
    fn test_normalize_frequencies() {
        let frequencies = [1, 2, 3, 0];
        assert_eq!(normalize_frequencies(&frequencies), [682, 1365, 2048, 0]);

        let frequencies = [0, 0, 0, 0];
        assert_eq!(normalize_frequencies(&frequencies), [0, 0, 0, 0]);
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
