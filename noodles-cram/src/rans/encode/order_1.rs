use std::{io, mem};

use byteorder::{LittleEndian, WriteBytesExt};

use super::{
    build_cumulative_frequencies, normalize, normalize_frequencies, update, BASE, LOWER_BOUND,
};

#[allow(dead_code)]
pub fn encode(data: &[u8]) -> io::Result<(Vec<Vec<u32>>, Vec<u8>)> {
    // Order-1 encoding does not support input smaller than 4 bytes.
    assert!(data.len() >= 4);

    let contexts = build_contexts(data, BASE);
    let freq = normalize_contexts(&contexts);
    let cfreq = build_cumulative_contexts(&freq);

    let mut buf = Vec::new();
    let mut states = [LOWER_BOUND; 4];

    // The input data is split into 4 equally sized chunks.
    let quarter = data.len() / states.len();

    let chunks = [
        &data[0..quarter],
        &data[quarter..2 * quarter],
        &data[2 * quarter..3 * quarter],
        &data[3 * quarter..4 * quarter],
    ];

    // The remainder of the input buffer is processed by the last state.
    if data.len() > 4 * quarter {
        // The last chunk includes the last symbol of the fourth chunk (the subtraction by 1). This
        // is safe because chunks are guaranteed to be nonempty.
        let remainder = &data[4 * quarter - 1..];

        for syms in remainder.windows(2).rev() {
            let freq_i = freq[syms[0] as usize][syms[1] as usize];
            let cfreq_i = cfreq[syms[0] as usize][syms[1] as usize];
            let x = normalize(&mut buf, states[3], freq_i)?;
            states[3] = update(x, freq_i, cfreq_i);
        }
    }

    for (((windows_0, windows_1), windows_2), windows_3) in chunks[0]
        .windows(2)
        .rev()
        .zip(chunks[1].windows(2).rev())
        .zip(chunks[2].windows(2).rev())
        .zip(chunks[3].windows(2).rev())
    {
        let windows = [windows_0, windows_1, windows_2, windows_3];

        for (state, syms) in states.iter_mut().zip(windows.iter()) {
            let freq_i = freq[syms[0] as usize][syms[1] as usize];
            let cfreq_i = cfreq[syms[0] as usize][syms[1] as usize];
            let x = normalize(&mut buf, *state, freq_i)?;
            *state = update(x, freq_i, cfreq_i);
        }
    }

    // The last state updates are for the starting contexts, i.e., `(0, chunks[i][0])`.
    for (state, chunk) in states.iter_mut().zip(chunks.iter()) {
        let sym = chunk[0] as usize;
        let freq_i = freq[0][sym];
        let cfreq_i = cfreq[0][sym];
        let x = normalize(&mut buf, *state, freq_i)?;
        *state = update(x, freq_i, cfreq_i);
    }

    let mut writer = Vec::with_capacity(4 * mem::size_of::<u32>() + buf.len());

    for &state in &states {
        writer.write_u32::<LittleEndian>(state)?;
    }

    writer.extend(buf.iter().rev());

    Ok((freq, writer))
}

fn build_contexts(data: &[u8], bin_count: usize) -> Vec<Vec<u32>> {
    let mut frequencies = vec![vec![0; bin_count]; bin_count];

    let quarter = data.len() / 4;

    for i in 0..4 {
        frequencies[0][data[i * quarter] as usize] += 1;
    }

    for window in data.windows(2) {
        let sym_0 = window[0] as usize;
        let sym_1 = window[1] as usize;
        frequencies[sym_0][sym_1] += 1;
    }

    frequencies
}

fn normalize_contexts(contexts: &[Vec<u32>]) -> Vec<Vec<u32>> {
    contexts
        .iter()
        .map(|frequencies| normalize_frequencies(frequencies))
        .collect()
}

fn build_cumulative_contexts(contexts: &[Vec<u32>]) -> Vec<Vec<u32>> {
    contexts
        .iter()
        .map(|frequencies| build_cumulative_frequencies(frequencies))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode() -> io::Result<()> {
        // § 14.4.2 Frequency table, Order-1 encoding (2020-07-22)
        let data = b"abracadabraabracadabraabracadabraabracadabra";
        let (actual_normalized_contexts, actual_compressed_data) = encode(data)?;

        let mut expected_normalized_contexts = [[0; BASE]; BASE];
        expected_normalized_contexts[0][usize::from(b'a')] = 4095;
        expected_normalized_contexts[usize::from(b'a')][usize::from(b'a')] = 646;
        expected_normalized_contexts[usize::from(b'a')][usize::from(b'b')] = 1725;
        expected_normalized_contexts[usize::from(b'a')][usize::from(b'c')] = 862;
        expected_normalized_contexts[usize::from(b'a')][usize::from(b'd')] = 862;
        expected_normalized_contexts[usize::from(b'b')][usize::from(b'r')] = 4095;
        expected_normalized_contexts[usize::from(b'c')][usize::from(b'a')] = 4095;
        expected_normalized_contexts[usize::from(b'd')][usize::from(b'a')] = 4095;
        expected_normalized_contexts[usize::from(b'r')][usize::from(b'a')] = 4095;

        assert_eq!(actual_normalized_contexts, expected_normalized_contexts);

        let expected_compressed_data = [
            0x75, 0x51, 0xc3, 0x3f, // states[0]
            0x75, 0x51, 0xc3, 0x3f, // states[1]
            0x75, 0x51, 0xc3, 0x3f, // states[2]
            0x75, 0x51, 0xc3, 0x3f, // states[3]
        ];

        assert_eq!(actual_compressed_data, expected_compressed_data);

        Ok(())
    }

    #[test]
    fn test_build_contexts() {
        let data = [1, 2, 3, 1, 2, 1, 2];
        let actual = build_contexts(&data, 4);
        let expected = [[0, 2, 1, 1], [0, 0, 3, 0], [0, 1, 0, 1], [0, 1, 0, 0]];
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_normalize_contexts() {
        let contexts = [
            vec![0, 2, 1, 1],
            vec![0, 0, 3, 0],
            vec![0, 1, 0, 1],
            vec![0, 1, 0, 0],
        ];
        let actual = normalize_contexts(&contexts);

        let expected = [
            [0, 2049, 1023, 1023],
            [0, 0, 4095, 0],
            [0, 2047, 0, 2048],
            [0, 4095, 0, 0],
        ];

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_build_cumulative_contexts() {
        let normalized_contexts = [
            vec![0, 2049, 1023, 1023],
            vec![0, 0, 4095, 0],
            vec![0, 2047, 0, 2048],
            vec![0, 4095, 0, 0],
        ];
        let actual = build_cumulative_contexts(&normalized_contexts);

        let expected = [
            [0, 0, 2049, 3072],
            [0, 0, 0, 4095],
            [0, 0, 2047, 2047],
            [0, 0, 4095, 4095],
        ];

        assert_eq!(actual, expected);
    }
}
