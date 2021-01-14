mod order_0;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use super::Order;

// Base `b`.
const BASE: usize = 256;

// Lower bound `L`.
const LOWER_BOUND: u32 = 0x800000;

#[allow(dead_code)]
pub fn rans_encode(order: Order, data: &[u8]) -> io::Result<Vec<u8>> {
    match order {
        Order::Zero => {
            let (normalized_frequencies, compressed_data) = order_0::encode(data)?;

            let mut writer = Vec::new();

            writer.write_u8(u8::from(Order::Zero))?;
            writer.write_u32::<LittleEndian>(compressed_data.len() as u32)?;
            writer.write_u32::<LittleEndian>(data.len() as u32)?;

            order_0::write_frequencies(&mut writer, &normalized_frequencies)?;

            writer.write_all(&compressed_data)?;

            Ok(writer)
        }
        Order::One => todo!(),
    }
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
}
