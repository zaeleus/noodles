use std::io::{self, Read};

use byteorder::ReadBytesExt;

use super::{order_0, read_states, state_cumulative_frequency, state_renormalize, state_step};
use crate::codecs::rans_4x8::{ALPHABET_SIZE, STATE_COUNT};

type Frequencies = [[u16; ALPHABET_SIZE]; ALPHABET_SIZE]; // F
type CumulativeFrequencies = Frequencies; // C
type CumulativeFrequenciesSymbolsTable = [[u8; 4096]; ALPHABET_SIZE];

pub fn decode<R>(reader: &mut R, dst: &mut [u8]) -> io::Result<()>
where
    R: Read,
{
    let frequencies = read_frequencies(reader)?;
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);

    let cumulative_frequencies_symbols_table =
        build_cumulative_frequencies_symbols_table(&cumulative_frequencies);

    let states = read_states(reader)?;

    let (chunk_0, chunk_1, chunk_2, chunk_3) = split_chunks(dst);
    let chunk_size = chunk_0.len();

    let mut chunks = [
        (states[0], 0, chunk_0),
        (states[1], 0, chunk_1),
        (states[2], 0, chunk_2),
        (states[3], 0, chunk_3),
    ];

    for k in 0..chunk_size {
        for (state, prev_sym, chunk) in &mut chunks {
            let i = usize::from(*prev_sym);
            let f = state_cumulative_frequency(*state);
            let sym = cumulative_frequencies_symbols_table[i][usize::from(f)];

            chunk[k] = sym;

            let j = usize::from(sym);
            *state = state_step(*state, frequencies[i][j], cumulative_frequencies[i][j]);
            *state = state_renormalize(*state, reader)?;

            *prev_sym = sym;
        }
    }

    let (mut state, mut prev_sym, chunk) = &mut chunks[3];
    let remainder = &mut chunk[chunk_size..];

    for d in remainder {
        let i = usize::from(prev_sym);
        let f = state_cumulative_frequency(state);
        let sym = cumulative_frequencies_symbols_table[i][usize::from(f)];

        *d = sym;

        let j = usize::from(sym);
        state = state_step(state, frequencies[i][j], cumulative_frequencies[i][j]);
        state = state_renormalize(state, reader)?;

        prev_sym = sym;
    }

    Ok(())
}

fn read_frequencies<R>(reader: &mut R) -> io::Result<Frequencies>
where
    R: Read,
{
    let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    let mut sym = reader.read_u8()?;
    let mut prev_sym = sym;

    loop {
        let f = order_0::read_frequencies(reader)?;
        frequencies[usize::from(sym)] = f;

        sym = reader.read_u8()?;

        if sym == 0 {
            break;
        }

        if sym - 1 == prev_sym {
            let len = reader.read_u8()?;

            for _ in 0..len {
                let f = order_0::read_frequencies(reader)?;
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

fn split_chunks(dst: &mut [u8]) -> (&mut [u8], &mut [u8], &mut [u8], &mut [u8]) {
    let chunk_size = dst.len() / STATE_COUNT;

    let (left_chunk, right_chunk) = dst.split_at_mut(2 * chunk_size);
    let (chunk_0, chunk_1) = left_chunk.split_at_mut(chunk_size);
    let (chunk_2, chunk_3) = right_chunk.split_at_mut(chunk_size);

    (chunk_0, chunk_1, chunk_2, chunk_3)
}
