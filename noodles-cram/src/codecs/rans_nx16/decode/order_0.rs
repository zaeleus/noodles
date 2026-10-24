use std::io;

use super::{SymbolTable, read_alphabet, read_states, renormalize_states, state_renormalize};
use crate::{codecs::rans_nx16::ALPHABET_SIZE, io::reader::num::read_uint7};

const NORMALIZATION_BITS: u32 = 12;

pub fn decode(src: &mut &[u8], dst: &mut [u8], state_count: usize) -> io::Result<()> {
    let frequencies = read_frequencies(src)?;
    // Small blocks retain the original lookup without allocating per-slot tables.
    if dst.len() < (1 << NORMALIZATION_BITS) {
        return match state_count {
            4 => decode_small::<4>(src, dst, &frequencies),
            32 => decode_small::<32>(src, dst, &frequencies),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid rANS Nx16 state count",
            )),
        };
    }
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

fn decode_small<const N: usize>(
    src: &mut &[u8],
    dst: &mut [u8],
    frequencies: &[u32; ALPHABET_SIZE],
) -> io::Result<()> {
    let mut cumulative = [0; ALPHABET_SIZE];
    let mut total = 0u64;
    for (c, &f) in cumulative.iter_mut().zip(frequencies) {
        *c = total as u32;
        total += u64::from(f);
        if total > (1 << NORMALIZATION_BITS) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid rANS Nx16 frequency table",
            ));
        }
    }
    let mut states = read_states::<N>(src)?;
    for chunk in dst.chunks_mut(N) {
        for (d, state) in chunk.iter_mut().zip(&mut states) {
            let f = *state & ((1 << NORMALIZATION_BITS) - 1);
            let mut sym = 0u8;
            while sym < u8::MAX && f >= cumulative[usize::from(sym + 1)] {
                sym += 1;
            }
            *d = sym;
            let i = usize::from(sym);
            *state = frequencies[i]
                .wrapping_mul(*state >> NORMALIZATION_BITS)
                .wrapping_add(f)
                .wrapping_sub(cumulative[i]);
            *state = state_renormalize(*state, src)?;
        }
    }
    Ok(())
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
