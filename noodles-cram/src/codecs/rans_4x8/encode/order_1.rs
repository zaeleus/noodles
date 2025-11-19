use std::io::{self, Write};

use super::{header, order_0, state_renormalize, state_step, write_header, write_states};
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

    let mut dst = vec![0; header::SIZE];

    write_frequencies(&mut dst, &frequencies)?;
    write_states(&mut dst, &states)?;
    dst.extend(buf.iter().rev());

    let compressed_size = dst[header::SIZE..].len();
    // SAFETY: `dst.len() >= header::SIZE`.
    let header_dst = dst.first_chunk_mut().unwrap();
    write_header(header_dst, Order::One, compressed_size, src.len())?;

    Ok(dst)
}

fn write_frequencies<W>(writer: &mut W, frequencies: &Frequencies) -> io::Result<()>
where
    W: Write,
{
    let alphabet = build_alphabet(frequencies);

    let mut iter = alphabet.iter().zip(frequencies).enumerate();
    let mut prev_sym = 0;

    while let Some((sym, (&a, f))) = iter.next() {
        if !a {
            continue;
        }

        // SAFETY: `sym <= ALPHABET_SIZE`.
        write_u8(writer, sym as u8)?;

        if sym > 0 && sym - 1 == prev_sym {
            let i = sym + 1;
            let len = alphabet[i..].iter().position(|&a| !a).unwrap_or(0);

            // SAFETY: `len < ALPHABET_SIZE`.
            write_u8(writer, len as u8)?;

            order_0::write_frequencies(writer, f)?;

            for (sym, (_, g)) in iter.by_ref().take(len) {
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

fn build_alphabet(frequencies: &Frequencies) -> [bool; ALPHABET_SIZE] {
    let mut alphabet = [false; ALPHABET_SIZE];

    for (a, fs) in alphabet.iter_mut().zip(frequencies) {
        *a = fs.iter().any(|&g| g > 0);
    }

    alphabet
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
    fn test_write_frequencies() -> io::Result<()> {
        // § 2.1.2 "Frequency table: Order-1 encoding" (563e8ab 2025-04-07):
        // "abracadabraabracadabraabracadabraabracadabrad".
        let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
        frequencies[usize::from(NUL)][usize::from(b'a')] = 4;
        frequencies[usize::from(b'a')][usize::from(b'a')] = 3;
        frequencies[usize::from(b'a')][usize::from(b'b')] = 8;
        frequencies[usize::from(b'a')][usize::from(b'c')] = 4;
        frequencies[usize::from(b'a')][usize::from(b'd')] = 5;
        frequencies[usize::from(b'b')][usize::from(b'r')] = 8;
        frequencies[usize::from(b'c')][usize::from(b'a')] = 4;
        frequencies[usize::from(b'd')][usize::from(b'a')] = 4;
        frequencies[usize::from(b'r')][usize::from(b'a')] = 8;

        let mut dst = Vec::new();
        write_frequencies(&mut dst, &frequencies)?;

        let expected = [
            0x00, // symbols[0] = '\0' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['\0']['a'] = 4
            0x00, // }
            b'a', // symbols[0] = 'a' {
            b'a', //   symbols[1] = 'a'
            0x03, //     frequencies['a']['a'] = 3
            b'b', //   symbols[1] = 'b'
            0x02, //     run length = 2
            0x08, //     frequencies['a']['b'] = 8
            0x04, //     frequencies['a']['c'] = 4
            0x05, //     frequencies['a']['d'] = 5
            0x00, // }
            b'b', // symbols[0] = 'b' {
            0x02, //   run length = 2
            b'r', //   symbols[1] = 'r'
            0x08, //     frequencies['b']['r'] = 8
            0x00, // }
            //    // symbols[0] = 'c' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['c']['a'] = 4
            0x00, // }
            //    // symbols[0] = 'd' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['d']['a'] = 4
            0x00, // }
            b'r', // symbols[0] = 'r' {
            b'a', //   symbols[1] = 'a'
            0x08, //     frequencies['r']['a'] = 8
            0x00, // }
            0x00, // EOF
        ];

        assert_eq!(dst, expected);

        Ok(())
    }

    #[test]
    fn test_build_alphabet() {
        // § 2.1.2 "Frequency table: Order-1 encoding" (563e8ab 2025-04-07)
        let src = b"abracadabraabracadabraabracadabraabracadabrad";
        let frequencies = build_raw_frequencies(src);
        let frequencies = normalize_frequencies(&frequencies);

        let mut expected = [false; ALPHABET_SIZE];
        expected[usize::from(NUL)] = true;
        expected[usize::from(b'a')] = true;
        expected[usize::from(b'b')] = true;
        expected[usize::from(b'c')] = true;
        expected[usize::from(b'd')] = true;
        expected[usize::from(b'r')] = true;

        assert_eq!(build_alphabet(&frequencies), expected);

        let src = b"noodles";
        let frequencies = build_raw_frequencies(src);
        let frequencies = normalize_frequencies(&frequencies);

        let mut expected = [false; ALPHABET_SIZE];
        expected[usize::from(NUL)] = true;
        expected[usize::from(b'n')] = true;
        expected[usize::from(b'o')] = true;
        expected[usize::from(b'd')] = true;
        expected[usize::from(b'l')] = true;
        expected[usize::from(b'e')] = true;

        assert_eq!(build_alphabet(&frequencies), expected);
    }

    #[test]
    fn test_build_raw_frequencies() {
        // § 2.1.2 "Frequency table: Order-1 encoding" (563e8ab 2025-04-07)
        let src = b"abracadabraabracadabraabracadabraabracadabrad";

        let mut expected = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
        expected[usize::from(NUL)][usize::from(b'a')] = 4;
        expected[usize::from(b'a')][usize::from(b'a')] = 3;
        expected[usize::from(b'a')][usize::from(b'b')] = 8;
        expected[usize::from(b'a')][usize::from(b'c')] = 4;
        expected[usize::from(b'a')][usize::from(b'd')] = 5;
        expected[usize::from(b'b')][usize::from(b'r')] = 8;
        expected[usize::from(b'c')][usize::from(b'a')] = 4;
        expected[usize::from(b'd')][usize::from(b'a')] = 4;
        expected[usize::from(b'r')][usize::from(b'a')] = 8;

        assert_eq!(build_raw_frequencies(src), expected);
    }

    #[test]
    fn test_normalize_frequencies() {
        // § 2.1.2 "Frequency table: Order-1 encoding" (563e8ab 2025-04-07)
        let src = b"abracadabraabracadabraabracadabraabracadabrad";
        let raw_frequencies = build_raw_frequencies(src);

        let mut expected = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
        expected[usize::from(NUL)][usize::from(b'a')] = 4095;
        expected[usize::from(b'a')][usize::from(b'a')] = 614;
        expected[usize::from(b'a')][usize::from(b'b')] = 1639;
        expected[usize::from(b'a')][usize::from(b'c')] = 819;
        expected[usize::from(b'a')][usize::from(b'd')] = 1023;
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

        // § 2.1.2 "Frequency table: Order-1 encoding" (563e8ab 2025-04-07)
        assert_eq!(
            split_chunks(b"abracadabraabracadabraabracadabraabracadabrad"),
            [
                &b"abracadabra"[..],
                &b"abracadabra"[..],
                &b"abracadabra"[..],
                &b"abracadabra"[..],
                &b"ad"[..],
            ]
        );
    }
}
