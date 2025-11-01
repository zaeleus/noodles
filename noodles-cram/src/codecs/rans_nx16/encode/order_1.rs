use std::io;

use super::{LOWER_BOUND, order_0, state_renormalize, state_step, write_alphabet, write_states};
use crate::{
    codecs::rans_nx16::ALPHABET_SIZE,
    io::writer::num::{write_u8, write_uint7},
};

const CONTEXT_SIZE: usize = 2;
const NORMALIZATION_BITS: u32 = 12;
const NUL: u8 = 0x00;

type Frequencies = [[u32; ALPHABET_SIZE]; ALPHABET_SIZE];

pub(super) struct Context {
    alphabet: [bool; ALPHABET_SIZE],
    frequencies: Frequencies,
}

pub(super) fn build_context(src: &[u8], state_count: usize) -> Context {
    let alphabet = build_alphabet(src);
    let frequencies = build_frequencies(src, state_count);
    let normalized_frequencies = normalize_frequencies(&frequencies);

    Context {
        alphabet,
        frequencies: normalized_frequencies,
    }
}

pub(super) fn write_context(dst: &mut Vec<u8>, ctx: &Context) -> io::Result<()> {
    // bits = 12, no compression (0)
    write_u8(dst, 12 << 4)?;
    write_alphabet(dst, &ctx.alphabet)?;
    write_frequencies(dst, &ctx.alphabet, &ctx.frequencies)?;
    Ok(())
}

pub fn encode(src: &[u8], ctx: &Context, state_count: usize, dst: &mut Vec<u8>) -> io::Result<()> {
    let frequencies = &ctx.frequencies;
    let cumulative_frequencies = build_cumulative_frequencies(frequencies);

    let mut buf = Vec::new();
    let mut states = vec![LOWER_BOUND; state_count];

    let (chunks, remainder) = split_chunks(src, state_count);

    if !remainder.is_empty() {
        let state = states.last_mut().unwrap();

        for syms in remainder.windows(CONTEXT_SIZE).rev() {
            let (i, j) = (usize::from(syms[0]), usize::from(syms[1]));
            let (f, g) = (frequencies[i][j], cumulative_frequencies[i][j]);
            *state = state_renormalize(*state, f, NORMALIZATION_BITS, &mut buf);
            *state = state_step(*state, f, g, NORMALIZATION_BITS);
        }
    }

    for (state, chunk) in states.iter_mut().rev().zip(chunks.iter().rev()) {
        for syms in chunk.windows(CONTEXT_SIZE).rev() {
            let (i, j) = (usize::from(syms[0]), usize::from(syms[1]));
            let (f, g) = (frequencies[i][j], cumulative_frequencies[i][j]);
            *state = state_renormalize(*state, f, NORMALIZATION_BITS, &mut buf);
            *state = state_step(*state, f, g, NORMALIZATION_BITS);
        }
    }

    for (state, chunk) in states.iter_mut().rev().zip(chunks.iter().rev()) {
        let (i, j) = (usize::from(NUL), usize::from(chunk[0]));
        let (f, g) = (frequencies[i][j], cumulative_frequencies[i][j]);
        *state = state_renormalize(*state, f, NORMALIZATION_BITS, &mut buf);
        *state = state_step(*state, f, g, NORMALIZATION_BITS);
    }

    write_states(dst, &states)?;
    dst.extend(buf.iter().rev());

    Ok(())
}

fn write_frequencies(
    dst: &mut Vec<u8>,
    alphabet: &[bool; ALPHABET_SIZE],
    frequencies: &Frequencies,
) -> io::Result<()> {
    for (_, fs) in alphabet.iter().zip(frequencies).filter(|(a, _)| **a) {
        let mut iter = alphabet.iter().zip(fs).filter(|(b, _)| **b).peekable();

        while let Some((_, &f)) = iter.next() {
            write_uint7(dst, f)?;

            if f == 0 {
                let mut len = 0;

                while let Some((_, g)) = iter.peek()
                    && **g == 0
                {
                    len += 1;
                    iter.next();
                }

                let n = u8::try_from(len)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                write_u8(dst, n)?;
            }
        }
    }

    Ok(())
}

fn build_alphabet(src: &[u8]) -> [bool; ALPHABET_SIZE] {
    let mut alphabet = [false; ALPHABET_SIZE];
    alphabet[usize::from(NUL)] = true;

    for &b in src {
        let i = usize::from(b);
        alphabet[i] = true;
    }

    alphabet
}

fn build_frequencies(src: &[u8], state_count: usize) -> [[u32; ALPHABET_SIZE]; ALPHABET_SIZE] {
    let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    let chunk_size = src.len() / state_count;

    for chunk in src.chunks_exact(chunk_size).take(state_count) {
        let (i, j) = (usize::from(NUL), usize::from(chunk[0]));
        frequencies[i][j] += 1;
    }

    for syms in src.windows(CONTEXT_SIZE) {
        let (i, j) = (usize::from(syms[0]), usize::from(syms[1]));
        frequencies[i][j] += 1;
    }

    frequencies
}

fn normalize_frequencies(raw_frequencies: &[[u32; ALPHABET_SIZE]; ALPHABET_SIZE]) -> Frequencies {
    let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    for (f, g) in raw_frequencies.iter().zip(&mut frequencies) {
        *g = order_0::normalize_frequencies(f);
    }

    frequencies
}

fn build_cumulative_frequencies(
    frequencies: &Frequencies,
) -> [[u32; ALPHABET_SIZE]; ALPHABET_SIZE] {
    let mut cumulative_frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    for (f, g) in frequencies.iter().zip(&mut cumulative_frequencies) {
        *g = order_0::build_cumulative_frequencies(f);
    }

    cumulative_frequencies
}

fn split_chunks(src: &[u8], state_count: usize) -> (Vec<&[u8]>, &[u8]) {
    let (q, r) = (src.len() / state_count, src.len() % state_count);
    let i = q * state_count;

    // SAFETY: `i` <= `src.len()`.
    let chunks = src[..i].chunks(q).collect();
    let remainder = if r > 0 { &src[i - 1..] } else { &[] };

    (chunks, remainder)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_frequencies() -> io::Result<()> {
        let src = b"abra";
        let frequencies = build_frequencies(src, 4);
        let alphabet = build_alphabet(src);

        let mut dst = Vec::new();
        write_frequencies(&mut dst, &alphabet, &frequencies)?;

        let expected = [
            0x00, // frequencies[NUL][NUL] = 0
            0x00, // run length = 0
            0x02, // frequencies[NUL]['a'] = 2
            0x01, // frequencies[NUL]['b'] = 1
            0x01, // frequencies[NUL]['r'] = 1
            0x00, // frequencies['a'][NUL] = 0
            0x01, // run length = 1
            0x01, // frequencies['a']['b'] = 1
            0x00, // frequencies['a']['r'] = 0
            0x00, // run length = 0
            0x00, // frequencies['b'][NUL] = 0
            0x02, // run length = 2
            0x01, // frequencies['b']['r'] = 1
            0x00, // frequencies['r'][NUL] = 0
            0x00, // run length = 0
            0x01, // frequencies['r']['a'] = 1
            0x00, // frequencies['r']['b'] = 0
            0x01, // run length = 1
        ];

        assert_eq!(dst, expected);

        Ok(())
    }

    #[test]
    fn test_build_alphabet() {
        let src = b"abracadabra";

        let mut expected = [false; ALPHABET_SIZE];
        expected[usize::from(NUL)] = true;
        expected[usize::from(b'a')] = true;
        expected[usize::from(b'b')] = true;
        expected[usize::from(b'c')] = true;
        expected[usize::from(b'd')] = true;
        expected[usize::from(b'r')] = true;

        assert_eq!(build_alphabet(src), expected);
    }

    #[test]
    fn test_build_frequencies() {
        let src = b"abracadabra";

        let mut expected = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
        expected[usize::from(NUL)][usize::from(b'a')] = 1;
        expected[usize::from(NUL)][usize::from(b'c')] = 1;
        expected[usize::from(NUL)][usize::from(b'd')] = 1;
        expected[usize::from(NUL)][usize::from(b'r')] = 1;
        expected[usize::from(b'a')][usize::from(b'b')] = 2;
        expected[usize::from(b'a')][usize::from(b'c')] = 1;
        expected[usize::from(b'a')][usize::from(b'd')] = 1;
        expected[usize::from(b'b')][usize::from(b'r')] = 2;
        expected[usize::from(b'c')][usize::from(b'a')] = 1;
        expected[usize::from(b'd')][usize::from(b'a')] = 1;
        expected[usize::from(b'r')][usize::from(b'a')] = 2;

        assert_eq!(build_frequencies(src, 4), expected);
    }

    #[test]
    fn test_split_chunks() {
        assert_eq!(
            split_chunks(b"abracadabra", 4),
            (
                vec![&b"ab"[..], &b"ra"[..], &b"ca"[..], &b"da"[..]],
                &b"abra"[..]
            )
        );

        assert_eq!(
            split_chunks(b"ndls", 2),
            (vec![&b"nd"[..], &b"ls"[..]], &[][..])
        );
    }
}
