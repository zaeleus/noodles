use std::io;

use super::{order_0, write_alphabet, write_states};
use crate::{
    codecs::rans_nx16::ALPHABET_SIZE,
    io::writer::num::{write_u8, write_uint7},
};

const CONTEXT_SIZE: usize = 2;

type Frequencies = [[u32; ALPHABET_SIZE]; ALPHABET_SIZE];

pub(super) struct Context {
    frequencies: Frequencies,
}

pub(super) fn build_context(src: &[u8], state_count: usize) -> Context {
    let frequencies = build_frequencies(src, state_count);
    let normalized_frequencies = normalize_frequencies(&frequencies);

    Context {
        frequencies: normalized_frequencies,
    }
}

pub(super) fn write_context(dst: &mut Vec<u8>, ctx: &Context) -> io::Result<()> {
    // bits = 12, no compression (0)
    write_u8(dst, 12 << 4)?;
    write_frequencies(dst, &ctx.frequencies)?;
    Ok(())
}

pub fn encode(src: &[u8], ctx: &Context, state_count: usize, dst: &mut Vec<u8>) -> io::Result<()> {
    use super::{LOWER_BOUND, normalize, update};

    let frequencies = &ctx.frequencies;
    let cumulative_frequencies = build_cumulative_frequencies(frequencies);

    let mut buf = Vec::new();
    let mut states = vec![LOWER_BOUND; state_count];

    let fraction = src.len() / state_count;

    if src.len() > state_count * fraction {
        let remainder = &src[state_count * fraction - 1..];

        for syms in remainder.windows(CONTEXT_SIZE).rev() {
            let (sym_0, sym_1) = (usize::from(syms[0]), usize::from(syms[1]));
            let freq_i = frequencies[sym_0][sym_1];
            let cfreq_i = cumulative_frequencies[sym_0][sym_1];
            let x = normalize(&mut buf, states[state_count - 1], freq_i, 12)?;
            states[state_count - 1] = update(x, cfreq_i, freq_i, 12);
        }
    }

    let chunks: Vec<_> = (0..state_count)
        .map(|i| &src[i * fraction..(i + 1) * fraction])
        .collect();

    let mut windows: Vec<_> = chunks
        .iter()
        .map(|chunk| chunk.windows(CONTEXT_SIZE).rev())
        .collect();

    let mut n = 0;
    let window_count = windows[0].size_hint().0;

    while n < window_count {
        for (state, ws) in states.iter_mut().rev().zip(windows.iter_mut().rev()) {
            let syms = ws.next().unwrap();
            let (sym_0, sym_1) = (usize::from(syms[0]), usize::from(syms[1]));
            let freq_i = frequencies[sym_0][sym_1];
            let cfreq_i = cumulative_frequencies[sym_0][sym_1];
            let x = normalize(&mut buf, *state, freq_i, 12)?;
            *state = update(x, cfreq_i, freq_i, 12);
        }

        n += 1;
    }

    for (state, chunk) in states.iter_mut().rev().zip(chunks.iter().rev()) {
        let sym = usize::from(chunk[0]);
        let freq_i = frequencies[0][sym];
        let cfreq_i = cumulative_frequencies[0][sym];
        let x = normalize(&mut buf, *state, freq_i, 12)?;
        *state = update(x, cfreq_i, freq_i, 12);
    }

    write_states(dst, &states)?;
    dst.extend(buf.iter().rev());

    Ok(())
}

fn write_frequencies(dst: &mut Vec<u8>, frequencies: &Frequencies) -> io::Result<()> {
    let mut alphabet = vec![0; 256];
    alphabet[0] = 1;

    for (i, f) in alphabet.iter_mut().enumerate() {
        for fs in frequencies {
            let g = fs[i];

            if g > 0 {
                *f += g;
            }
        }
    }

    write_alphabet(dst, &alphabet)?;

    for (sym_0, fs) in frequencies.iter().enumerate() {
        if alphabet[sym_0] == 0 {
            continue;
        }

        let mut rle = 0;

        for (sym_1, &f) in fs.iter().enumerate() {
            if alphabet[sym_1] == 0 {
                continue;
            }

            if rle > 0 {
                rle -= 1;
            } else {
                write_uint7(dst, f)?;

                if f == 0 {
                    for (sym, &g) in alphabet.iter().enumerate().skip(sym_1 + 1) {
                        if g == 0 {
                            continue;
                        }

                        if frequencies[sym_0][sym] == 0 {
                            rle += 1;
                        } else {
                            break;
                        }
                    }

                    write_u8(dst, rle)?;
                }
            }
        }
    }

    Ok(())
}

fn build_frequencies(src: &[u8], state_count: usize) -> [[u32; ALPHABET_SIZE]; ALPHABET_SIZE] {
    let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    let fraction = src.len() / state_count;

    for i in 0..state_count {
        let sym = usize::from(src[i * fraction]);
        frequencies[0][sym] += 1;
    }

    for window in src.windows(CONTEXT_SIZE) {
        let sym_0 = usize::from(window[0]);
        let sym_1 = usize::from(window[1]);
        frequencies[sym_0][sym_1] += 1;
    }

    let sym = src.last().copied().map(usize::from).unwrap();
    frequencies[sym][0] += 1;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_frequencies() {
        const NUL: u8 = 0x00;

        let src = b"abracadabra";

        let mut expected = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
        expected[usize::from(NUL)][usize::from(b'a')] = 1;
        expected[usize::from(NUL)][usize::from(b'c')] = 1;
        expected[usize::from(NUL)][usize::from(b'd')] = 1;
        expected[usize::from(NUL)][usize::from(b'r')] = 1;
        expected[usize::from(b'a')][usize::from(NUL)] = 1;
        expected[usize::from(b'a')][usize::from(b'b')] = 2;
        expected[usize::from(b'a')][usize::from(b'c')] = 1;
        expected[usize::from(b'a')][usize::from(b'd')] = 1;
        expected[usize::from(b'b')][usize::from(b'r')] = 2;
        expected[usize::from(b'c')][usize::from(b'a')] = 1;
        expected[usize::from(b'd')][usize::from(b'a')] = 1;
        expected[usize::from(b'r')][usize::from(b'a')] = 2;

        assert_eq!(build_frequencies(src, 4), expected);
    }
}
