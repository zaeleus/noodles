mod order_0;
mod order_1;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use super::Order;

// Base `b`.
const BASE: usize = 256;

// Lower bound `L`.
const LOWER_BOUND: u32 = 0x800000;

#[allow(dead_code)]
pub fn rans_encode(order: Order, data: &[u8]) -> io::Result<Vec<u8>> {
    let compressed_blob = match order {
        Order::Zero => {
            let (normalized_frequencies, compressed_data) = order_0::encode(data)?;

            let mut compressed_blob = Vec::new();
            order_0::write_frequencies(&mut compressed_blob, &normalized_frequencies)?;
            compressed_blob.extend(&compressed_data);

            compressed_blob
        }
        Order::One => {
            let (normalized_contexts, compressed_data) = order_1::encode(data)?;

            let mut compressed_blob = Vec::new();
            order_1::write_contexts(&mut compressed_blob, &normalized_contexts)?;
            compressed_blob.extend(&compressed_data);

            compressed_blob
        }
    };

    let mut writer = Vec::new();

    write_header(
        &mut writer,
        order,
        compressed_blob.len() as u32,
        data.len() as u32,
    )?;

    writer.write_all(&compressed_blob)?;

    Ok(writer)
}

fn write_header<W>(
    writer: &mut W,
    order: Order,
    compressed_len: u32,
    data_len: u32,
) -> io::Result<()>
where
    W: Write,
{
    writer.write_u8(u8::from(order))?;
    writer.write_u32::<LittleEndian>(compressed_len)?;
    writer.write_u32::<LittleEndian>(data_len)?;
    Ok(())
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

    if sum == 0 {
        return vec![0; frequencies.len()];
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
    fn test_rans_encode_with_order_0() -> io::Result<()> {
        let data = b"noodles";
        let actual = rans_encode(Order::Zero, data)?;

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
    fn test_rans_encode_with_order_1() -> io::Result<()> {
        let data = b"noodles";
        let actual = rans_encode(Order::One, data)?;

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
