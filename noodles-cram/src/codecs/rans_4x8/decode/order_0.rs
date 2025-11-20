use std::io;

use super::{read_states, state_cumulative_frequency, state_renormalize, state_step};
use crate::{
    codecs::rans_4x8::ALPHABET_SIZE,
    io::reader::num::{read_itf8_as, read_u8},
};

type Frequencies = [u16; ALPHABET_SIZE]; // F
type CumulativeFrequencies = Frequencies; // C
type CumulativeFrequenciesSymbolsTable = [u8; 4096];

pub fn decode(src: &mut &[u8], dst: &mut [u8]) -> io::Result<()> {
    let frequencies = read_frequencies(src)?;
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);

    let cumulative_frequencies_symbols_table =
        build_cumulative_frequencies_symbols_table(&cumulative_frequencies);

    let mut states = read_states(src)?;

    for chunk in dst.chunks_mut(states.len()) {
        for (d, state) in chunk.iter_mut().zip(states.iter_mut()) {
            let f = state_cumulative_frequency(*state);
            let sym = cumulative_frequencies_symbols_table[usize::from(f)];

            *d = sym;

            let i = usize::from(sym);
            *state = state_step(*state, frequencies[i], cumulative_frequencies[i]);
            *state = state_renormalize(*state, src)?;
        }
    }

    Ok(())
}

pub(super) fn read_frequencies(src: &mut &[u8]) -> io::Result<Frequencies> {
    const NUL: u8 = 0x00;

    let mut frequencies = [0; ALPHABET_SIZE];

    let mut sym = read_u8(src)?;
    let mut prev_sym = sym;

    loop {
        let f = read_itf8_as(src)?;
        frequencies[usize::from(sym)] = f;

        sym = read_u8(src)?;

        if sym == NUL {
            break;
        }

        if sym - 1 == prev_sym {
            let len = read_u8(src)?;

            for _ in 0..len {
                let f = read_itf8_as(src)?;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_frequencies() -> io::Result<()> {
        let src = [
            b'a', // symbol = 'a'
            0x05, // frequencies['a'] = 5
            b'b', // symbol = 'b'
            0x02, // run length = 2
            0x02, // frequencies['b'] = 2
            0x01, // frequencies['c'] = 1
            0x01, // frequencies['d'] = 1
            b'r', // symbol = 'r'
            0x02, // frequencies['r'] = 2
            0x00, // EOF
        ];

        let mut expected = [0; ALPHABET_SIZE];
        expected[usize::from(b'a')] = 5;
        expected[usize::from(b'b')] = 2;
        expected[usize::from(b'c')] = 1;
        expected[usize::from(b'd')] = 1;
        expected[usize::from(b'r')] = 2;

        assert_eq!(read_frequencies(&mut &src[..])?, expected);

        Ok(())
    }

    #[test]
    fn test_build_cumulative_frequencies() {
        let mut frequencies = [0; ALPHABET_SIZE];
        frequencies[usize::from(b'a')] = 5;
        frequencies[usize::from(b'b')] = 2;
        frequencies[usize::from(b'c')] = 1;
        frequencies[usize::from(b'd')] = 1;
        frequencies[usize::from(b'r')] = 2;

        let mut expected = [0; ALPHABET_SIZE];
        expected[..usize::from(b'b')].fill(0);
        expected[usize::from(b'b')..usize::from(b'c')].fill(5);
        expected[usize::from(b'c')..usize::from(b'd')].fill(7);
        expected[usize::from(b'd')..usize::from(b'e')].fill(8);
        expected[usize::from(b'e')..usize::from(b's')].fill(9);
        expected[usize::from(b's')..].fill(11);

        assert_eq!(build_cumulative_frequencies(&frequencies), expected);
    }
}
