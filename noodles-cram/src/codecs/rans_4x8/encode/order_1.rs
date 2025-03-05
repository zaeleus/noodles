use std::io::{self, Write};

use byteorder::WriteBytesExt;

use super::{normalize, order_0, update, write_states};
use crate::codecs::rans_4x8::{ALPHABET_SIZE, LOWER_BOUND, STATE_COUNT};

const NUL: u8 = 0x00;

pub fn encode(src: &[u8]) -> io::Result<Vec<u8>> {
    use super::{write_header, Order};

    // Order-1 encoding does not support input smaller than 4 bytes.
    assert!(src.len() >= 4);

    let raw_frequencies = build_raw_frequencies(src);
    let frequencies = normalize_frequencies(&raw_frequencies);
    let cfreq = build_cumulative_contexts(&frequencies);

    let mut buf = Vec::new();
    let mut states = [LOWER_BOUND; STATE_COUNT];

    // The input data is split into 4 equally sized chunks.
    let quarter = src.len() / states.len();

    let chunks = [
        &src[0..quarter],
        &src[quarter..2 * quarter],
        &src[2 * quarter..3 * quarter],
        &src[3 * quarter..4 * quarter],
    ];

    // The remainder of the input buffer is processed by the last state.
    if src.len() > 4 * quarter {
        // The last chunk includes the last symbol of the fourth chunk (the subtraction by 1). This
        // is safe because chunks are guaranteed to be nonempty.
        let remainder = &src[4 * quarter - 1..];

        for syms in remainder.windows(2).rev() {
            let freq_i = frequencies[usize::from(syms[0])][usize::from(syms[1])];
            let cfreq_i = cfreq[usize::from(syms[0])][usize::from(syms[1])];
            let x = normalize(&mut buf, states[3], freq_i)?;
            states[3] = update(x, freq_i, cfreq_i);
        }
    }

    let mut windows = [
        chunks[0].windows(2).rev(),
        chunks[1].windows(2).rev(),
        chunks[2].windows(2).rev(),
        chunks[3].windows(2).rev(),
    ];

    let mut n = 0;
    let window_count = windows[0].size_hint().0;

    while n < window_count {
        for (state, ws) in states.iter_mut().rev().zip(windows.iter_mut().rev()) {
            let syms = ws.next().unwrap();
            let freq_i = frequencies[usize::from(syms[0])][usize::from(syms[1])];
            let cfreq_i = cfreq[usize::from(syms[0])][usize::from(syms[1])];
            let x = normalize(&mut buf, *state, freq_i)?;
            *state = update(x, freq_i, cfreq_i);
        }

        n += 1;
    }

    // The last state updates are for the starting contexts, i.e., `(0, chunks[i][0])`.
    for (state, chunk) in states.iter_mut().rev().zip(chunks.iter().rev()) {
        let sym = usize::from(chunk[0]);
        let freq_i = frequencies[0][sym];
        let cfreq_i = cfreq[0][sym];
        let x = normalize(&mut buf, *state, freq_i)?;
        *state = update(x, freq_i, cfreq_i);
    }

    let mut dst = vec![0; 9];

    write_contexts(&mut dst, &frequencies)?;
    write_states(&mut dst, &states)?;
    dst.extend(buf.iter().rev());

    let compressed_len = dst[9..].len();
    let mut writer = &mut dst[..9];
    write_header(&mut writer, Order::One, compressed_len, src.len())?;

    Ok(dst)
}

fn write_contexts<W>(
    writer: &mut W,
    contexts: &[[u16; ALPHABET_SIZE]; ALPHABET_SIZE],
) -> io::Result<()>
where
    W: Write,
{
    use super::order_0;

    let sums: Vec<u16> = contexts
        .iter()
        .map(|frequencies| frequencies.iter().sum())
        .collect();

    let mut rle = 0;

    for (sym, (frequencies, &sum)) in contexts.iter().zip(sums.iter()).enumerate() {
        if sum == 0 {
            continue;
        }

        if rle > 0 {
            rle -= 1;
        } else {
            writer.write_u8(sym as u8)?;

            if sym > 0 && sums[sym - 1] > 0 {
                rle = sums[sym + 1..].iter().position(|&s| s == 0).unwrap_or(0);
                writer.write_u8(rle as u8)?;
            }
        }

        order_0::write_frequencies(writer, frequencies)?;
    }

    writer.write_u8(0x00)?;

    Ok(())
}

fn build_raw_frequencies(src: &[u8]) -> [[u32; ALPHABET_SIZE]; ALPHABET_SIZE] {
    const CONTEXT_SIZE: usize = 2;

    assert!(src.len() >= STATE_COUNT);

    let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    // § 2.1.2 "Frequency table: Order-1 encoding": "We use the ASCII value (`\0`) as the starting
    // context for each interleaved rANS state..."
    let chunk_size = src.len() / STATE_COUNT;

    for chunk in src.chunks_exact(chunk_size).take(STATE_COUNT) {
        // SAFETY: `chunk.len() > 0`.
        let (i, j) = (usize::from(NUL), usize::from(chunk[0]));
        frequencies[i][j] += 1;
    }

    for syms in src.windows(CONTEXT_SIZE) {
        let (i, j) = (usize::from(syms[0]), usize::from(syms[1]));
        frequencies[i][j] += 1;
    }

    frequencies
}

fn normalize_frequencies(
    raw_frequencies: &[[u32; ALPHABET_SIZE]; ALPHABET_SIZE],
) -> [[u16; ALPHABET_SIZE]; ALPHABET_SIZE] {
    let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    for (f, g) in raw_frequencies.iter().zip(&mut frequencies) {
        *g = order_0::normalize_frequencies(f);
    }

    frequencies
}

fn build_cumulative_contexts(contexts: &[[u16; ALPHABET_SIZE]; ALPHABET_SIZE]) -> Vec<Vec<u16>> {
    contexts
        .iter()
        .map(|frequencies| order_0::build_cumulative_frequencies(frequencies).to_vec())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode() -> io::Result<()> {
        // § 14.4.2 Frequency table, Order-1 encoding (2020-07-22)
        let data = b"abracadabraabracadabraabracadabraabracadabra";
        let actual = encode(data)?;

        let expected = [
            // header
            0x01, // order = 1
            0x36, 0x00, 0x00, 0x00, // compressed_len = 54
            0x2c, 0x00, 0x00, 0x00, // uncompressed_len = 44
            // frequency table
            0x00, // syms[0] = '\0' {
            0x61, // syms[1] = 'a'
            0x8f, 0xff, // f[0]['a'] = 4095
            0x00, // }
            0x61, // syms[0] = 'a' {
            0x61, // syms[1] = 'a'
            0x82, 0x86, // f['a']['a'] = 646
            0x62, // syms[1] = 'b'
            0x02, // rle = 2
            0x86, 0xbd, // f['a']['b'] = 1725
            0x83, 0x5e, // f['a']['c'] = 862
            0x83, 0x5e, // f['a']['d'] = 862
            0x00, // }
            0x62, // syms[0] = 'b' {
            0x02, // rle = 2
            0x72, // syms[1] = 'r'
            0x8f, 0xff, // f['b']['r'] = 4095
            0x00, // }
            // syms[0] = 'c' {
            0x61, // 'a'
            0x8f, 0xff, // f['c']['a'] = 4095
            0x00, // }
            // syms[0] = 'd' {
            0x61, // syms[1] = 'a'
            0x8f, 0xff, // f['d']['a'] = 4095
            0x00, // }
            0x72, // syms[0] = 'r' {
            0x61, // syms[1] = 'a'
            0x8f, 0xff, // f['r']['a'] = 4095
            0x00, // }
            0x00, // end
            // compressed blob
            0x75, 0x51, 0xc3, 0x3f, // states[0]
            0x75, 0x51, 0xc3, 0x3f, // states[1]
            0x75, 0x51, 0xc3, 0x3f, // states[2]
            0x75, 0x51, 0xc3, 0x3f, // states[3]
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_build_raw_frequencies() {
        // § 2.1.2 "Frequency table: Order-1 encoding" (2023-03-15)
        let src = b"abracadabraabracadabraabracadabraabracadabr";

        let mut expected = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
        expected[usize::from(NUL)][usize::from(b'a')] = 2;
        expected[usize::from(NUL)][usize::from(b'r')] = 1;
        expected[usize::from(NUL)][usize::from(b'b')] = 1;
        expected[usize::from(b'a')][usize::from(b'a')] = 3;
        expected[usize::from(b'a')][usize::from(b'b')] = 8;
        expected[usize::from(b'a')][usize::from(b'c')] = 4;
        expected[usize::from(b'a')][usize::from(b'd')] = 4;
        expected[usize::from(b'b')][usize::from(b'r')] = 8;
        expected[usize::from(b'c')][usize::from(b'a')] = 4;
        expected[usize::from(b'd')][usize::from(b'a')] = 4;
        expected[usize::from(b'r')][usize::from(b'a')] = 7;

        assert_eq!(build_raw_frequencies(src), expected);
    }

    #[test]
    fn test_normalize_contexts() {
        // § 2.1.2 "Frequency table: Order-1 encoding" (2023-03-15)
        // let src = b"abracadabraabracadabraabracadabraabracadabr";

        let src = b"abracadabraabracadabraabracadabraabracadabra";
        let raw_frequencies = build_raw_frequencies(src);

        let mut expected = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
        expected[usize::from(NUL)][usize::from(b'a')] = 4095;
        expected[usize::from(b'a')][usize::from(b'a')] = 646;
        expected[usize::from(b'a')][usize::from(b'b')] = 1725;
        expected[usize::from(b'a')][usize::from(b'c')] = 862;
        expected[usize::from(b'a')][usize::from(b'd')] = 862;
        expected[usize::from(b'b')][usize::from(b'r')] = 4095;
        expected[usize::from(b'c')][usize::from(b'a')] = 4095;
        expected[usize::from(b'd')][usize::from(b'a')] = 4095;
        expected[usize::from(b'r')][usize::from(b'a')] = 4095;

        assert_eq!(normalize_frequencies(&raw_frequencies), expected);
    }
}
