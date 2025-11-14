use std::io;

use super::{order_0, read_states, state_cumulative_frequency, state_renormalize, state_step};
use crate::{
    codecs::rans_4x8::{ALPHABET_SIZE, STATE_COUNT},
    io::reader::num::read_u8,
};

type Frequencies = [[u16; ALPHABET_SIZE]; ALPHABET_SIZE]; // F
type CumulativeFrequencies = Frequencies; // C
type CumulativeFrequenciesSymbolsTable = [[u8; 4096]; ALPHABET_SIZE];

pub fn decode(src: &mut &[u8], dst: &mut [u8]) -> io::Result<()> {
    let frequencies = read_frequencies(src)?;
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);

    let cumulative_frequencies_symbols_table =
        build_cumulative_frequencies_symbols_table(&cumulative_frequencies);

    let mut states = read_states(src)?;
    let mut prev_syms = [0; STATE_COUNT];

    let [chunk_0, chunk_1, chunk_2, chunk_3, chunk_4] = split_chunks(dst);
    let chunks = chunk_0.iter_mut().zip(chunk_1).zip(chunk_2).zip(chunk_3);

    for (((d0, d1), d2), d3) in chunks {
        let dsts = [d0, d1, d2, d3];

        for (state, (prev_sym, d)) in states.iter_mut().zip(prev_syms.iter_mut().zip(dsts)) {
            let i = usize::from(*prev_sym);
            let f = state_cumulative_frequency(*state);
            let sym = cumulative_frequencies_symbols_table[i][usize::from(f)];

            *d = sym;

            let j = usize::from(sym);
            *state = state_step(*state, frequencies[i][j], cumulative_frequencies[i][j]);
            *state = state_renormalize(*state, src)?;

            *prev_sym = sym;
        }
    }

    let mut state = states[3];
    let mut prev_sym = prev_syms[3];

    for d in chunk_4 {
        let i = usize::from(prev_sym);
        let f = state_cumulative_frequency(state);
        let sym = cumulative_frequencies_symbols_table[i][usize::from(f)];

        *d = sym;

        let j = usize::from(sym);
        state = state_step(state, frequencies[i][j], cumulative_frequencies[i][j]);
        state = state_renormalize(state, src)?;

        prev_sym = sym;
    }

    Ok(())
}

fn read_frequencies(src: &mut &[u8]) -> io::Result<Frequencies> {
    let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    let mut sym = read_u8(src)?;
    let mut prev_sym = sym;

    loop {
        let f = order_0::read_frequencies(src)?;
        frequencies[usize::from(sym)] = f;

        sym = read_u8(src)?;

        if sym == 0 {
            break;
        }

        if sym - 1 == prev_sym {
            let len = read_u8(src)?;

            for _ in 0..len {
                let f = order_0::read_frequencies(src)?;
                frequencies[usize::from(sym)] = f;
                sym += 1;
            }
        }

        prev_sym = sym;
    }

    Ok(frequencies)
}

fn build_cumulative_frequencies(frequencies: &Frequencies) -> CumulativeFrequencies {
    let mut cumulative_frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    for (f, g) in frequencies.iter().zip(&mut cumulative_frequencies) {
        *g = order_0::build_cumulative_frequencies(f);
    }

    cumulative_frequencies
}

pub fn build_cumulative_frequencies_symbols_table(
    cumulative_freqs: &CumulativeFrequencies,
) -> Box<CumulativeFrequenciesSymbolsTable> {
    let mut tables = Box::new([[0; 4096]; 256]);

    for (table, cumulative_freqs) in tables.iter_mut().zip(cumulative_freqs) {
        *table = order_0::build_cumulative_frequencies_symbols_table(cumulative_freqs);
    }

    tables
}

fn split_chunks(dst: &mut [u8]) -> [&mut [u8]; 5] {
    let chunk_size = dst.len() / STATE_COUNT;

    let (left_chunk, right_chunk) = dst.split_at_mut(2 * chunk_size);
    let (chunk_0, chunk_1) = left_chunk.split_at_mut(chunk_size);
    let (chunk_2, chunk_3_4) = right_chunk.split_at_mut(chunk_size);
    let (chunk_3, chunk_4) = chunk_3_4.split_at_mut(chunk_size);

    [chunk_0, chunk_1, chunk_2, chunk_3, chunk_4]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_frequencies() -> io::Result<()> {
        const NUL: u8 = 0x00;

        // § 2.1.2 "Frequency table: Order-1 encoding" (563e8ab 2025-04-07):
        // "abracadabraabracadabraabracadabraabracadabrad".
        let src = [
            0x00, // symbols[0] = '\0' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['\0']['a'] = 4
            0x00, // }
            b'a', // symbols[0] = 'a' {
            b'a', //   symbols[1] = 'a'
            0x03, //     frequencies['a']['a'] = 3
            b'b', //   symbols[1] = 'b'
            0x02, //     run length = 2
            0x08, //     frequencies['a']['b'] = 8
            0x04, //     frequencies['a']['c'] = 4
            0x05, //     frequencies['a']['d'] = 5
            0x00, // }
            b'b', // symbols[0] = 'b' {
            0x02, //   run length = 2
            b'r', //   symbols[1] = 'r'
            0x08, //     frequencies['b']['r'] = 8
            0x00, // }
            //    // symbols[0] = 'c' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['c']['a'] = 4
            0x00, // }
            //    // symbols[0] = 'd' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['d']['a'] = 4
            0x00, // }
            b'r', // symbols[0] = 'r' {
            b'a', //   symbols[1] = 'a'
            0x08, //     frequencies['r']['a'] = 8
            0x00, // }
            0x00, // EOF
        ];

        let mut expected = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
        expected[usize::from(NUL)][usize::from(b'a')] = 4;
        expected[usize::from(b'a')][usize::from(b'a')] = 3;
        expected[usize::from(b'a')][usize::from(b'b')] = 8;
        expected[usize::from(b'a')][usize::from(b'c')] = 4;
        expected[usize::from(b'a')][usize::from(b'd')] = 5;
        expected[usize::from(b'b')][usize::from(b'r')] = 8;
        expected[usize::from(b'c')][usize::from(b'a')] = 4;
        expected[usize::from(b'd')][usize::from(b'a')] = 4;
        expected[usize::from(b'r')][usize::from(b'a')] = 8;

        assert_eq!(read_frequencies(&mut &src[..])?, expected);

        Ok(())
    }
}
