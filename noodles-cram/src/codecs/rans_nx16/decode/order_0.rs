use std::io;

use super::{SymbolTable, read_alphabet, read_states, renormalize_states, state_renormalize};
use crate::{codecs::rans_nx16::ALPHABET_SIZE, io::reader::num::read_uint7};

const NORMALIZATION_BITS: u32 = 12;

pub fn decode(src: &mut &[u8], dst: &mut [u8], state_count: usize) -> io::Result<()> {
    let frequencies = read_frequencies(src)?;
    let table = SymbolTable::new(&frequencies, NORMALIZATION_BITS)?;

    match state_count {
        4 => decode_with_states::<4>(src, dst, &table),
        32 => decode_with_states::<32>(src, dst, &table),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid rANS Nx16 state count",
        )),
    }
}

fn decode_with_states<const N: usize>(
    src: &mut &[u8],
    dst: &mut [u8],
    table: &SymbolTable,
) -> io::Result<()> {
    let mut states = read_states::<N>(src)?;

    let (chunks, tail) = dst.as_chunks_mut::<N>();

    for chunk in chunks {
        for (d, state) in chunk.iter_mut().zip(&mut states) {
            *d = table.decode_symbol(state, NORMALIZATION_BITS);
        }

        // The state steps are independent, so they run ahead of the ordered renormalization reads.
        renormalize_states(&mut states, src)?;
    }

    for (d, state) in tail.iter_mut().zip(&mut states) {
        *d = table.decode_symbol(state, NORMALIZATION_BITS);
        *state = state_renormalize(*state, src)?;
    }

    Ok(())
}

pub fn normalize_frequencies(frequencies: &mut [u32; ALPHABET_SIZE], bits: u32) {
    // Sum in u64: a malformed table must not overflow before the symbol table rejects it.
    let mut sum: u64 = frequencies.iter().copied().map(u64::from).sum();

    if sum == 0 || sum >= (1 << bits) {
        return;
    }

    let mut shift = 0;

    while sum < (1 << bits) {
        sum *= 2;
        shift += 1;
    }

    for f in frequencies {
        *f <<= shift;
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
