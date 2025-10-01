use std::io;

use super::{
    cumulative_frequencies_symbol, order_0, read_alphabet, read_states, split_off,
    state_cumulative_frequency, state_renormalize, state_step,
};
use crate::{
    codecs::rans_nx16::ALPHABET_SIZE,
    io::reader::num::{read_u8, read_uint7, read_uint7_as},
};

type Frequencies = [[u32; ALPHABET_SIZE]; ALPHABET_SIZE];
type CumulativeFrequencies = Frequencies;

pub fn decode(src: &mut &[u8], dst: &mut [u8], state_count: usize) -> io::Result<()> {
    let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
    let bits = read_frequencies(src, &mut frequencies)?;

    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);

    let mut states = read_states(src, state_count)?;

    let mut i = 0;
    let mut last_syms = vec![0; states.len()];

    while i < dst.len() / state_count {
        for j in 0..state_count {
            let f = state_cumulative_frequency(states[j], bits);
            let s = cumulative_frequencies_symbol(&cumulative_frequencies[last_syms[j]], f);

            dst[i + j * (dst.len() / state_count)] = s;

            states[j] = state_step(
                states[j],
                frequencies[last_syms[j]][usize::from(s)],
                cumulative_frequencies[last_syms[j]][usize::from(s)],
                bits,
            );

            states[j] = state_renormalize(states[j], src)?;

            last_syms[j] = usize::from(s);
        }

        i += 1;
    }

    i *= state_count;
    let m = state_count - 1;

    while i < dst.len() {
        let f = state_cumulative_frequency(states[m], bits);
        let s = cumulative_frequencies_symbol(&cumulative_frequencies[last_syms[m]], f);

        dst[i] = s;

        states[m] = state_step(
            states[m],
            frequencies[last_syms[m]][usize::from(s)],
            cumulative_frequencies[last_syms[m]][usize::from(s)],
            bits,
        );

        states[m] = state_renormalize(states[m], src)?;

        last_syms[m] = usize::from(s);

        i += 1;
    }

    Ok(())
}

fn read_frequencies(src: &mut &[u8], frequencies: &mut Frequencies) -> io::Result<u32> {
    const STATE_COUNT: usize = 4;

    let n = read_u8(src)?;
    let bits = u32::from(n >> 4);
    let is_compressed = (n & 0x01) != 0;

    if is_compressed {
        let uncompressed_size = read_uint7_as(src)?;
        let compressed_size = read_uint7_as(src)?;
        let mut compressed_data = split_off(src, compressed_size)?;
        let mut dst = vec![0; uncompressed_size];
        order_0::decode(&mut compressed_data, &mut dst, STATE_COUNT)?;
        read_frequencies_inner(&mut &dst[..], frequencies, bits)?;
    } else {
        read_frequencies_inner(src, frequencies, bits)?;
    }

    Ok(bits)
}

fn read_frequencies_inner(
    src: &mut &[u8],
    frequencies: &mut Frequencies,
    bits: u32,
) -> io::Result<()> {
    let alphabet = read_alphabet(src)?;

    for (_, fs) in alphabet.iter().zip(frequencies).filter(|(a, _)| **a) {
        let mut iter = alphabet.iter().zip(fs.iter_mut()).filter(|(b, _)| **b);

        while let Some((_, f)) = iter.next() {
            *f = read_uint7(src)?;

            if *f == 0 {
                let n = read_u8(src).map(usize::from)?;
                for _ in iter.by_ref().take(n) {}
            }
        }

        order_0::normalize_frequencies(fs, bits);
    }

    Ok(())
}

fn build_cumulative_frequencies(frequencies: &Frequencies) -> CumulativeFrequencies {
    let mut cumulative_frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    for (f, g) in frequencies.iter().zip(&mut cumulative_frequencies) {
        *g = order_0::build_cumulative_frequencies(f);
    }

    cumulative_frequencies
}
