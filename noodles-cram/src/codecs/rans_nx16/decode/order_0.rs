use std::io;

use super::{
    cumulative_frequencies_symbol, read_alphabet, read_states, state_cumulative_frequency,
    state_renormalize, state_step,
};
use crate::{codecs::rans_nx16::ALPHABET_SIZE, io::reader::num::read_uint7};

const NORMALIZATION_BITS: u32 = 12;

pub fn decode(src: &mut &[u8], dst: &mut [u8], state_count: usize) -> io::Result<()> {
    let frequencies = read_frequencies(src)?;
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);

    let mut states = read_states(src, state_count)?;

    for (i, d) in dst.iter_mut().enumerate() {
        let j = i % state_count;

        let f = state_cumulative_frequency(states[j], NORMALIZATION_BITS);
        let s = cumulative_frequencies_symbol(&cumulative_frequencies, f);

        *d = s;

        states[j] = state_step(
            states[j],
            frequencies[s as usize],
            cumulative_frequencies[s as usize],
            NORMALIZATION_BITS,
        );

        states[j] = state_renormalize(states[j], src)?;
    }

    Ok(())
}

pub fn normalize_frequencies(freqs: &mut [u32], bits: u32) {
    let mut total: u32 = freqs.iter().sum();

    if total == 0 || total == (1 << bits) {
        return;
    }

    let mut shift = 0;

    while total < (1 << bits) {
        total *= 2;
        shift += 1;
    }

    for freq in freqs {
        *freq <<= shift;
    }
}

fn read_frequencies(src: &mut &[u8]) -> io::Result<[u32; ALPHABET_SIZE]> {
    let alphabet = read_alphabet(src)?;

    let mut frequencies = [0; ALPHABET_SIZE];

    for (i, frequency) in alphabet.iter().zip(&mut frequencies) {
        if *i {
            *frequency = read_uint7(src)?;
        }
    }

    normalize_frequencies(&mut frequencies, NORMALIZATION_BITS);

    Ok(frequencies)
}

fn build_cumulative_frequencies(frequencies: &[u32; ALPHABET_SIZE]) -> [u32; ALPHABET_SIZE] {
    let mut cumulative_frequencies = [0; ALPHABET_SIZE];

    let mut f = cumulative_frequencies[0];

    for (next_f, g) in cumulative_frequencies[1..].iter_mut().zip(frequencies) {
        *next_f = f + g;
        f = *next_f;
    }

    cumulative_frequencies
}
