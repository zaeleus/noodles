use std::io::{self, Read};

use super::{
    cumulative_frequencies_symbol, order_0, read_states, state_cumulative_frequency,
    state_renormalize, state_step,
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
    let comp = read_u8(src)?;
    let bits = u32::from(comp >> 4);

    if comp & 0x01 != 0 {
        let u_size = read_uint7_as(src)?;
        let c_size = read_uint7_as(src)?;

        let mut c_data = vec![0; c_size];
        src.read_exact(&mut c_data)?;

        let mut c_data_reader = &c_data[..];
        let mut u_data = vec![0; u_size];
        order_0::decode(&mut c_data_reader, &mut u_data, 4)?;

        let mut u_data_reader = &u_data[..];
        read_frequencies_inner(&mut u_data_reader, frequencies, bits)?;
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
    use super::read_alphabet;

    let alphabet = read_alphabet(src)?;

    for (i, a) in alphabet.iter().enumerate() {
        if !a {
            continue;
        }

        let mut run = 0;

        for (j, b) in alphabet.iter().enumerate() {
            if !b {
                continue;
            }

            if run > 0 {
                run -= 1;
            } else {
                let f = read_uint7(src)?;

                frequencies[i][j] = f;

                if f == 0 {
                    run = read_u8(src)?;
                }
            }
        }

        order_0::normalize_frequencies(&mut frequencies[i], bits);
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
