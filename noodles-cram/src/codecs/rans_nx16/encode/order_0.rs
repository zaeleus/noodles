use std::{
    io::{self, Write},
    mem,
};

use super::{build_frequencies, normalize_frequencies};
use crate::{
    codecs::rans_nx16::ALPHABET_SIZE,
    io::writer::num::{write_u32_le, write_uint7},
};

const NORMALIZATION_BITS: u32 = 12;

pub(super) struct Context {
    pub(super) frequencies: [u32; ALPHABET_SIZE],
}

pub(super) fn build_context(src: &[u8]) -> Context {
    let frequencies = build_frequencies(src);
    let normalized_frequencies = normalize_frequencies(&frequencies);

    Context {
        frequencies: normalized_frequencies,
    }
}

pub fn encode(src: &[u8], ctx: &Context, n: usize) -> io::Result<Vec<u8>> {
    use super::{LOWER_BOUND, build_cumulative_frequencies, normalize, update};

    let frequencies = &ctx.frequencies;
    let cumulative_frequencies = build_cumulative_frequencies(frequencies);

    let mut buf = Vec::new();
    let mut states = vec![LOWER_BOUND; n];

    for (i, &sym) in src.iter().enumerate().rev() {
        let j = i % states.len();

        let mut x = states[j];
        let freq_i = frequencies[usize::from(sym)];
        let cfreq_i = cumulative_frequencies[usize::from(sym)];

        x = normalize(&mut buf, x, freq_i, NORMALIZATION_BITS)?;
        states[j] = update(x, cfreq_i, freq_i, NORMALIZATION_BITS);
    }

    let mut dst = Vec::with_capacity(states.len() * mem::size_of::<u32>() + buf.len());

    for &state in &states {
        write_u32_le(&mut dst, state)?;
    }

    dst.extend(buf.iter().rev());

    Ok(dst)
}

pub fn write_frequencies<W>(writer: &mut W, frequencies: &[u32]) -> io::Result<()>
where
    W: Write,
{
    use super::write_alphabet;

    write_alphabet(writer, frequencies)?;

    for &f in frequencies {
        if f > 0 {
            write_uint7(writer, f)?;
        }
    }

    Ok(())
}
