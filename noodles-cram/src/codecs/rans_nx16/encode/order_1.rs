use std::{
    io::{self, Write},
    mem,
};

use super::{order_0, write_states};
use crate::io::writer::num::{write_u8, write_uint7};

pub fn encode(src: &[u8], n: usize) -> io::Result<(Vec<Vec<u32>>, Vec<u8>)> {
    use super::{LOWER_BOUND, normalize, update};

    let contexts = build_contexts(src, n);
    let freq = normalize_contexts(contexts);
    let cfreq = build_cumulative_contexts(&freq);

    let mut buf = Vec::new();
    let mut states = vec![LOWER_BOUND; n];

    let fraction = src.len() / n;

    if src.len() > n * fraction {
        let remainder = &src[n * fraction - 1..];

        for syms in remainder.windows(2).rev() {
            let (sym_0, sym_1) = (usize::from(syms[0]), usize::from(syms[1]));
            let freq_i = freq[sym_0][sym_1];
            let cfreq_i = cfreq[sym_0][sym_1];
            let x = normalize(&mut buf, states[n - 1], freq_i, 12)?;
            states[n - 1] = update(x, cfreq_i, freq_i, 12);
        }
    }

    let chunks: Vec<_> = (0..n)
        .map(|i| &src[i * fraction..(i + 1) * fraction])
        .collect();

    let mut windows: Vec<_> = chunks.iter().map(|chunk| chunk.windows(2).rev()).collect();

    let mut n = 0;
    let window_count = windows[0].size_hint().0;

    while n < window_count {
        for (state, ws) in states.iter_mut().rev().zip(windows.iter_mut().rev()) {
            let syms = ws.next().unwrap();
            let (sym_0, sym_1) = (usize::from(syms[0]), usize::from(syms[1]));
            let freq_i = freq[sym_0][sym_1];
            let cfreq_i = cfreq[sym_0][sym_1];
            let x = normalize(&mut buf, *state, freq_i, 12)?;
            *state = update(x, cfreq_i, freq_i, 12);
        }

        n += 1;
    }

    for (state, chunk) in states.iter_mut().rev().zip(chunks.iter().rev()) {
        let sym = usize::from(chunk[0]);
        let freq_i = freq[0][sym];
        let cfreq_i = cfreq[0][sym];
        let x = normalize(&mut buf, *state, freq_i, 12)?;
        *state = update(x, cfreq_i, freq_i, 12);
    }

    let mut dst = Vec::with_capacity(n * mem::size_of::<u32>() + buf.len());
    write_states(&mut dst, &states)?;
    dst.extend(buf.iter().rev());

    Ok((freq, dst))
}

pub fn write_contexts<W>(writer: &mut W, contexts: &[Vec<u32>]) -> io::Result<()>
where
    W: Write,
{
    use super::write_alphabet;

    let mut frequencies = vec![0; 256];
    frequencies[0] = 1;

    for (i, f) in frequencies.iter_mut().enumerate() {
        for context in contexts {
            let g = context[i];

            if g > 0 {
                *f += g;
            }
        }
    }

    write_alphabet(writer, &frequencies)?;

    for (sym_0, context) in contexts.iter().enumerate() {
        if frequencies[sym_0] == 0 {
            continue;
        }

        let mut rle = 0;

        for (sym_1, &f) in context.iter().enumerate() {
            if frequencies[sym_1] == 0 {
                continue;
            }

            if rle > 0 {
                rle -= 1;
            } else {
                write_uint7(writer, f)?;

                if f == 0 {
                    for (sym, &g) in frequencies.iter().enumerate().skip(sym_1 + 1) {
                        if g == 0 {
                            continue;
                        }

                        if contexts[sym_0][sym] == 0 {
                            rle += 1;
                        } else {
                            break;
                        }
                    }

                    write_u8(writer, rle)?;
                }
            }
        }
    }

    Ok(())
}

fn build_contexts(src: &[u8], n: usize) -> Vec<Vec<u32>> {
    let mut frequencies = vec![vec![0; 256]; 256];

    let fraction = src.len() / n;

    for i in 0..n {
        let sym = usize::from(src[i * fraction]);
        frequencies[0][sym] += 1;
    }

    for window in src.windows(2) {
        let sym_0 = usize::from(window[0]);
        let sym_1 = usize::from(window[1]);
        frequencies[sym_0][sym_1] += 1;
    }

    let sym = src.last().copied().map(usize::from).unwrap();
    frequencies[sym][0] += 1;

    frequencies
}

fn normalize_contexts(contexts: Vec<Vec<u32>>) -> Vec<Vec<u32>> {
    contexts
        .into_iter()
        .map(|frequencies| order_0::normalize_frequencies(&frequencies).to_vec())
        .collect()
}

fn build_cumulative_contexts(contexts: &[Vec<u32>]) -> Vec<Vec<u32>> {
    contexts
        .iter()
        .map(|frequencies| order_0::build_cumulative_frequencies(frequencies))
        .collect()
}
