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
    let mut prev_syms = vec![0; state_count];
    let mut i = 0;

    while i < dst.len() / state_count {
        for (j, (state, prev_sym)) in states.iter_mut().zip(&mut prev_syms).enumerate() {
            let k = usize::from(*prev_sym);

            let f = state_cumulative_frequency(*state, bits);
            let sym = cumulative_frequencies_symbol(&cumulative_frequencies[k], f);

            dst[j * (dst.len() / state_count) + i] = sym;

            let l = usize::from(sym);

            *state = state_step(
                *state,
                frequencies[k][l],
                cumulative_frequencies[k][l],
                bits,
            );

            *state = state_renormalize(*state, src)?;

            *prev_sym = sym;
        }

        i += 1;
    }

    i *= state_count;

    let m = state_count - 1;
    let mut state = states[m];
    let mut prev_sym = prev_syms[m];

    while i < dst.len() {
        let k = usize::from(prev_sym);
        let f = state_cumulative_frequency(state, bits);
        let sym = cumulative_frequencies_symbol(&cumulative_frequencies[k], f);

        dst[i] = sym;

        let l = usize::from(sym);
        state = state_step(state, frequencies[k][l], cumulative_frequencies[k][l], bits);
        state = state_renormalize(state, src)?;

        prev_sym = sym;

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
