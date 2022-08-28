use std::{
    io::{self, Write},
    mem,
};

use byteorder::{LittleEndian, WriteBytesExt};

use crate::writer::num::write_uint7;

pub fn encode(src: &[u8], n: usize) -> io::Result<(Vec<Vec<u32>>, Vec<u8>)> {
    use super::{normalize, update};

    const LOWER_BOUND: u32 = 0x8000;

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

    let windows: Vec<_> = chunks.iter().map(|chunk| chunk.windows(2).rev()).collect();

    for ws in windows {
        for (state, syms) in states.iter_mut().zip(ws) {
            let (sym_0, sym_1) = (usize::from(syms[0]), usize::from(syms[1]));
            let freq_i = freq[sym_0][sym_1];
            let cfreq_i = cfreq[sym_0][sym_1];
            let x = normalize(&mut buf, *state, freq_i, 12)?;
            *state = update(x, cfreq_i, freq_i, 12);
        }
    }

    for (state, chunk) in states.iter_mut().zip(chunks.iter()) {
        let sym = usize::from(chunk[0]);
        let freq_i = freq[0][sym];
        let cfreq_i = cfreq[0][sym];
        let x = normalize(&mut buf, *state, freq_i, 12)?;
        *state = update(x, cfreq_i, freq_i, 12);
    }

    let mut dst = Vec::with_capacity(n * mem::size_of::<u32>() + buf.len());

    for &state in &states {
        dst.write_u32::<LittleEndian>(state)?;
    }

    dst.extend(buf.iter().rev());

    Ok((freq, dst))
}

pub fn write_contexts<W>(writer: &mut W, contexts: &[Vec<u32>]) -> io::Result<()>
where
    W: Write,
{
    use super::write_alphabet;

    let sums: Vec<u32> = contexts
        .iter()
        .map(|frequencies| frequencies.iter().sum())
        .collect();

    write_alphabet(writer, &sums)?;

    for (sym_0, context) in contexts.iter().enumerate() {
        if sums[sym_0] == 0 {
            continue;
        }

        let mut rle = 0;

        for (sym_1, &f) in context.iter().enumerate() {
            if sums[sym_1] == 0 {
                continue;
            }

            if rle > 0 {
                rle -= 1;
            } else {
                write_uint7(writer, f)?;

                if f == 0 {
                    for (sym, &sum) in sums.iter().enumerate().skip(sym_1 + 1) {
                        if sum == 0 {
                            continue;
                        }

                        if contexts[sym_0][sym] == 0 {
                            rle += 1;
                        } else {
                            break;
                        }
                    }

                    writer.write_u8(rle)?;
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

    frequencies
}

fn normalize_contexts(contexts: Vec<Vec<u32>>) -> Vec<Vec<u32>> {
    use super::normalize_frequencies;

    contexts
        .into_iter()
        .map(|frequencies| normalize_frequencies(&frequencies))
        .collect()
}

fn build_cumulative_contexts(contexts: &[Vec<u32>]) -> Vec<Vec<u32>> {
    use super::build_cumulative_frequencies;

    contexts
        .iter()
        .map(|frequencies| build_cumulative_frequencies(frequencies))
        .collect()
}
