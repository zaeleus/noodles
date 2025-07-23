use std::io::{self, Read};

use super::{read_states, state_cumulative_frequency, state_renormalize, state_step};
use crate::{
    codecs::rans_4x8::ALPHABET_SIZE,
    io::reader::num::{read_itf8_as, read_u8},
};

type Frequencies = [u16; ALPHABET_SIZE]; // F
type CumulativeFrequencies = Frequencies; // C
type CumulativeFrequenciesSymbolsTable = [u8; 4096];

pub fn decode<R>(reader: &mut R, dst: &mut [u8]) -> io::Result<()>
where
    R: Read,
{
    let frequencies = read_frequencies(reader)?;
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);

    let cumulative_frequencies_symbols_table =
        build_cumulative_frequencies_symbols_table(&cumulative_frequencies);

    let mut states = read_states(reader)?;

    for chunk in dst.chunks_mut(states.len()) {
        for (d, state) in chunk.iter_mut().zip(states.iter_mut()) {
            let f = state_cumulative_frequency(*state);
            let sym = cumulative_frequencies_symbols_table[usize::from(f)];

            *d = sym;

            let i = usize::from(sym);
            *state = state_step(*state, frequencies[i], cumulative_frequencies[i]);
            *state = state_renormalize(*state, reader)?;
        }
    }

    Ok(())
}

pub(super) fn read_frequencies<R>(reader: &mut R) -> io::Result<Frequencies>
where
    R: Read,
{
    const NUL: u8 = 0x00;

    let mut frequencies = [0; ALPHABET_SIZE];

    let mut sym = read_u8(reader)?;
    let mut prev_sym = sym;

    loop {
        let f = read_itf8_as(reader)?;
        frequencies[usize::from(sym)] = f;

        sym = read_u8(reader)?;

        if sym == NUL {
            break;
        }

        if sym - 1 == prev_sym {
            let len = read_u8(reader)?;

            for _ in 0..len {
                let f = read_itf8_as(reader)?;
                frequencies[usize::from(sym)] = f;
                sym += 1;
            }
        }

        prev_sym = sym;
    }

    Ok(frequencies)
}

pub(super) fn build_cumulative_frequencies(frequencies: &Frequencies) -> CumulativeFrequencies {
    let mut cumulative_frequencies = [0; ALPHABET_SIZE];

    let mut f = cumulative_frequencies[0];

    for (next_f, g) in cumulative_frequencies[1..].iter_mut().zip(frequencies) {
        *next_f = f + g;
        f = *next_f;
    }

    cumulative_frequencies
}

pub(super) fn build_cumulative_frequencies_symbols_table(
    cumulative_freqs: &CumulativeFrequencies,
) -> CumulativeFrequenciesSymbolsTable {
    let mut table = [0; 4096];
    let mut sym = 0;

    for (f, g) in (0u16..).zip(&mut table) {
        while sym < u8::MAX && f >= cumulative_freqs[usize::from(sym + 1)] {
            sym += 1;
        }

        *g = sym;
    }

    table
}
