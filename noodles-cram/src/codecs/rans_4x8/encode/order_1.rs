use std::io::{self, Write};

use super::{order_0, state_renormalize, state_step, write_header, write_states};
use crate::{
    codecs::rans_4x8::{ALPHABET_SIZE, LOWER_BOUND, Order, STATE_COUNT},
    io::writer::num::write_u8,
};

const CONTEXT_SIZE: usize = 2;
const NUL: u8 = 0x00;

type RawFrequencies = [[u32; ALPHABET_SIZE]; ALPHABET_SIZE];
type Frequencies = [[u16; ALPHABET_SIZE]; ALPHABET_SIZE]; // F
type CumulativeFrequencies = Frequencies; // C

pub fn encode(src: &[u8]) -> io::Result<Vec<u8>> {
    // § 2.2.1 "rANS entropy encoding: Interleaving" (2023-03-15): "We do not permit Order-1
    // encoding of data streams smaller than 4 bytes."
    if src.len() < STATE_COUNT {
        return Err(io::Error::from(io::ErrorKind::InvalidInput));
    }

    let raw_frequencies = build_raw_frequencies(src);
    let frequencies = normalize_frequencies(&raw_frequencies);
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);

    let mut buf = Vec::new();
    let mut states = [LOWER_BOUND; STATE_COUNT];

    let [chunk_0, chunk_1, chunk_2, chunk_3, chunk_4] = split_chunks(src);

    // § 2.2.1 "rANS entropy encoding: Interleaving" (2023-03-15): "Any remainder, when the input
    // buffer is not divisible by 4, is processed ... by the 4th rANS state."
    for syms in chunk_4.windows(CONTEXT_SIZE).rev() {
        let (i, j) = (usize::from(syms[0]), usize::from(syms[1]));
        states[3] = state_renormalize(states[3], frequencies[i][j], &mut buf)?;
        states[3] = state_step(states[3], frequencies[i][j], cumulative_frequencies[i][j]);
    }

    let mut windows_rev = chunk_0
        .windows(CONTEXT_SIZE)
        .rev()
        .zip(chunk_1.windows(CONTEXT_SIZE).rev())
        .zip(chunk_2.windows(CONTEXT_SIZE).rev())
        .zip(chunk_3.windows(CONTEXT_SIZE).rev());

    for (((w0, w1), w2), w3) in &mut windows_rev {
        let windows = [w0, w1, w2, w3];

        for (state, syms) in states.iter_mut().rev().zip(windows.iter().rev()) {
            let (i, j) = (usize::from(syms[0]), usize::from(syms[1]));
            *state = state_renormalize(*state, frequencies[i][j], &mut buf)?;
            *state = state_step(*state, frequencies[i][j], cumulative_frequencies[i][j]);
        }
    }

    let chunks = [chunk_0, chunk_1, chunk_2, chunk_3];

    for (state, chunk) in states.iter_mut().rev().zip(chunks.iter().rev()) {
        let (i, j) = (usize::from(NUL), usize::from(chunk[0]));
        *state = state_renormalize(*state, frequencies[i][j], &mut buf)?;
        *state = state_step(*state, frequencies[i][j], cumulative_frequencies[i][j]);
    }

    let mut dst = vec![0; 9];

    write_frequencies(&mut dst, &frequencies)?;
    write_states(&mut dst, &states)?;
    dst.extend(buf.iter().rev());

    let compressed_size = dst[9..].len();
    let mut writer = &mut dst[..9];
    write_header(&mut writer, Order::One, compressed_size, src.len())?;

    Ok(dst)
}

fn write_frequencies<W>(writer: &mut W, frequencies: &Frequencies) -> io::Result<()>
where
    W: Write,
{
    let mut statuses = [false; ALPHABET_SIZE];

    for (s, f) in statuses.iter_mut().zip(frequencies) {
        *s = !f.iter().any(|&g| g > 0);
    }

    let mut iter = frequencies.iter().zip(&statuses).enumerate();
    let mut prev_sym = 0;

    while let Some((sym, (f, &is_empty))) = iter.next() {
        if is_empty {
            continue;
        }

        // SAFETY: `sym <= ALPHABET_SIZE`.
        write_u8(writer, sym as u8)?;

        if sym > 0 && sym - 1 == prev_sym {
            let i = sym + 1;
            let len = statuses[i..].iter().position(|s| *s).unwrap_or(0);

            // SAFETY: `len < ALPHABET_SIZE`.
            write_u8(writer, len as u8)?;

            order_0::write_frequencies(writer, f)?;

            for (sym, (g, _)) in iter.by_ref().take(len) {
                order_0::write_frequencies(writer, g)?;
                prev_sym = sym;
            }

            continue;
        }

        order_0::write_frequencies(writer, f)?;

        prev_sym = sym;
    }

    write_u8(writer, NUL)?;

    Ok(())
}

fn build_raw_frequencies(src: &[u8]) -> RawFrequencies {
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

fn normalize_frequencies(raw_frequencies: &RawFrequencies) -> Frequencies {
    let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    for (f, g) in raw_frequencies.iter().zip(&mut frequencies) {
        *g = order_0::normalize_frequencies(f);
    }

    frequencies
}

fn build_cumulative_frequencies(frequencies: &Frequencies) -> CumulativeFrequencies {
    let mut cumulative_frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    for (f, g) in frequencies.iter().zip(&mut cumulative_frequencies) {
        *g = order_0::build_cumulative_frequencies(f);
    }

    cumulative_frequencies
}

fn split_chunks(buf: &[u8]) -> [&[u8]; 5] {
    let chunk_size = buf.len() / STATE_COUNT;

    let (left_chunk, right_chunk) = buf.split_at(2 * chunk_size);
    let (chunk_0, chunk_1) = left_chunk.split_at(chunk_size);
    let (chunk_2, chunk_3_4) = right_chunk.split_at(chunk_size);
    let (chunk_3, _) = chunk_3_4.split_at(chunk_size);

    // The last chunk includes the last byte of chunk 3 for context.
    let chunk_4 = &chunk_3_4[chunk_size - 1..];

    [chunk_0, chunk_1, chunk_2, chunk_3, chunk_4]
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
            0x36, 0x00, 0x00, 0x00, // compressed size = 54
            0x2c, 0x00, 0x00, 0x00, // uncompressed size = 44
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
    fn test_normalize_frequencies() {
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

    #[test]
    fn test_split_chunks() {
        assert_eq!(
            split_chunks(&[1, 2, 3, 4]),
            [&[1][..], &[2][..], &[3][..], &[4][..], &[4][..]]
        );

        assert_eq!(
            split_chunks(&[1, 2, 3, 4, 5]),
            [&[1][..], &[2][..], &[3][..], &[4][..], &[4, 5][..]]
        );
    }
}
