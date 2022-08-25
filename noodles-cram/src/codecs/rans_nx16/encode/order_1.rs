use std::io::{self, Write};

use byteorder::WriteBytesExt;

use crate::writer::num::write_uint7;

pub fn encode(src: &[u8], n: usize) -> io::Result<(Vec<Vec<u32>>, Vec<u8>)> {
    let contexts = build_contexts(src, n);
    let freq = normalize_contexts(contexts);
    let _cfreq = build_cumulative_contexts(&freq);

    Ok((freq, Vec::new()))
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

    let quarter = src.len() / n;

    for i in 0..n {
        let sym = usize::from(src[i * quarter]);
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
