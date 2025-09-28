use std::io;

use super::{read_alphabet, read_states};
use crate::io::reader::num::read_uint7;

pub fn decode(src: &mut &[u8], dst: &mut [u8], state_count: usize) -> io::Result<()> {
    use super::{
        rans_advance_step_nx16, rans_get_cumulative_freq_nx16, rans_get_symbol_from_freq,
        rans_renorm_nx16,
    };

    let frequencies = read_frequencies(src)?;
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);

    let mut states = read_states(src, state_count)?;

    for (i, d) in dst.iter_mut().enumerate() {
        let j = i % state_count;

        let f = rans_get_cumulative_freq_nx16(states[j], 12);
        let s = rans_get_symbol_from_freq(&cumulative_frequencies, f);

        *d = s;

        states[j] = rans_advance_step_nx16(
            states[j],
            cumulative_frequencies[s as usize],
            frequencies[s as usize],
            12,
        );

        states[j] = rans_renorm_nx16(src, states[j])?;
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

fn read_frequencies(src: &mut &[u8]) -> io::Result<[u32; 256]> {
    let mut frequencies = [0; 256];

    let alphabet = read_alphabet(src)?;

    for i in 0..alphabet.len() {
        if alphabet[i] {
            frequencies[i] = read_uint7(src)?;
        }
    }

    normalize_frequencies(&mut frequencies, 12);

    Ok(frequencies)
}

fn build_cumulative_frequencies(frequencies: &[u32; 256]) -> [u32; 256] {
    let mut cumulative_frequencies = [0; 256];

    for i in 0..255 {
        cumulative_frequencies[i + 1] = cumulative_frequencies[i] + frequencies[i];
    }

    cumulative_frequencies
}
