use std::io;

use super::{
    apply_step, build_symbol_table, order_0, read_alphabet, read_states, renormalize_states,
    split_off, state_renormalize,
};
use crate::{
    codecs::rans_nx16::ALPHABET_SIZE,
    io::reader::num::{read_u8, read_uint7, read_uint7_as},
};

/// The largest frequency table shift.
///
/// htscodecs writes 10 or 12. Every context with a nonzero row gets a `1 << bits` table, so a
/// larger shift would only let a few bytes of malformed input allocate and fill tens of MiB.
const MAX_BITS: u32 = 12;

/// Symbol tables (see [`SymbolTable`]) for every context, stored contiguously so that a symbol
/// costs one context offset lookup rather than a pointer chase. Context `c` has the `1 << bits`
/// slots starting at `offsets[c]` in both `symbols` and `steps`. Contexts absent from the
/// alphabet, and contexts whose frequencies are all zero, share an empty table at offset 0.
struct Tables {
    symbols: Vec<u8>,
    steps: Vec<u32>,
    offsets: [usize; ALPHABET_SIZE],
}

impl Tables {
    #[inline(always)]
    fn decode_symbol(&self, prev_sym: u8, state: &mut u32, bits: u32) -> u8 {
        let i = self.offsets[usize::from(prev_sym)] + (*state & ((1 << bits) - 1)) as usize;
        let sym = self.symbols[i];
        *state = apply_step(*state, self.steps[i], bits);
        sym
    }
}

pub fn decode(src: &mut &[u8], dst: &mut [u8], state_count: usize) -> io::Result<()> {
    let (tables, bits) = read_frequencies(src)?;

    match state_count {
        4 => decode_with_states::<4>(src, dst, &tables, bits),
        32 => decode_with_states::<32>(src, dst, &tables, bits),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid rANS Nx16 state count",
        )),
    }
}

fn decode_with_states<const N: usize>(
    src: &mut &[u8],
    dst: &mut [u8],
    tables: &Tables,
    bits: u32,
) -> io::Result<()> {
    let mut states = read_states::<N>(src)?;
    let mut prev_syms = [0u8; N];

    let chunk_size = dst.len() / N;
    let (body, last_chunk) = dst.split_at_mut(chunk_size * N);

    if chunk_size > 0 {
        let mut chunks = body.chunks_exact_mut(chunk_size);
        // `body` has exactly `N` chunks of `chunk_size`.
        let mut dsts: [&mut [u8]; N] = std::array::from_fn(|_| chunks.next().unwrap());

        for i in 0..chunk_size {
            for ((d, state), prev_sym) in dsts.iter_mut().zip(&mut states).zip(&mut prev_syms) {
                let sym = tables.decode_symbol(*prev_sym, state, bits);
                d[i] = sym;
                *prev_sym = sym;
            }

            // The state steps are independent, so they run ahead of the ordered renormalization
            // reads.
            renormalize_states(&mut states, src)?;
        }
    }

    let mut state = states[N - 1];
    let mut prev_sym = prev_syms[N - 1];

    for d in last_chunk {
        let sym = tables.decode_symbol(prev_sym, &mut state, bits);
        *d = sym;
        state = state_renormalize(state, src)?;
        prev_sym = sym;
    }

    Ok(())
}

fn read_frequencies(src: &mut &[u8]) -> io::Result<(Tables, u32)> {
    const STATE_COUNT: usize = 4;

    let n = read_u8(src)?;
    let bits = u32::from(n >> 4);
    let is_compressed = (n & 0x01) != 0;

    if bits > MAX_BITS {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid rANS Nx16 frequency table",
        ));
    }

    let tables = if is_compressed {
        let uncompressed_size = read_uint7_as(src)?;
        let compressed_size = read_uint7_as(src)?;
        let mut compressed_data = split_off(src, compressed_size)?;
        let mut dst = vec![0; uncompressed_size];
        order_0::decode(&mut compressed_data, &mut dst, STATE_COUNT)?;
        read_frequencies_inner(&mut &dst[..], bits)?
    } else {
        read_frequencies_inner(src, bits)?
    };

    Ok((tables, bits))
}

fn read_frequencies_inner(src: &mut &[u8], bits: u32) -> io::Result<Tables> {
    let alphabet = read_alphabet(src)?;

    // Every row is read before any table is allocated, so all-zero rows take no space.
    let mut rows = Vec::new();

    for (context, _) in alphabet.iter().enumerate().filter(|(_, a)| **a) {
        let mut fs = [0; ALPHABET_SIZE];
        let mut iter = alphabet.iter().zip(fs.iter_mut()).filter(|(b, _)| **b);

        while let Some((_, f)) = iter.next() {
            *f = read_uint7(src)?;

            if *f == 0 {
                let n = read_u8(src).map(usize::from)?;
                for _ in iter.by_ref().take(n) {}
            }
        }

        if fs.iter().any(|&f| f > 0) {
            order_0::normalize_frequencies(&mut fs, bits);
            rows.push((context, fs));
        }
    }

    let size = 1 << bits;

    let mut symbols = vec![0; (rows.len() + 1) * size];
    let mut steps = vec![0; (rows.len() + 1) * size];
    let mut offsets = [0; ALPHABET_SIZE];

    build_symbol_table(
        &[0; ALPHABET_SIZE],
        bits,
        &mut symbols[..size],
        &mut steps[..size],
    )?;

    for (i, (context, fs)) in rows.iter().enumerate() {
        let offset = (i + 1) * size;
        let range = offset..offset + size;
        build_symbol_table(fs, bits, &mut symbols[range.clone()], &mut steps[range])?;
        offsets[*context] = offset;
    }

    Ok(Tables {
        symbols,
        steps,
        offsets,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_frequencies_shares_the_empty_table_and_rejects_shifts_above_12() -> io::Result<()>
    {
        // Shift 12, uncompressed; alphabet {A, C}; row A: A=1, C=1; row C: all zero (0, then a
        // run of 1 more zero).
        let mut src = [0xc0, b'A', b'C', 0x00, 1, 1, 0, 1];
        let (tables, bits) = read_frequencies(&mut &src[..])?;

        assert_eq!(bits, 12);
        assert_eq!(tables.symbols.len(), 2 << 12);
        assert_eq!(tables.offsets[usize::from(b'A')], 1 << 12);
        assert_eq!(tables.offsets[usize::from(b'C')], 0);

        src[0] = 0xd0; // shift 13
        let result = read_frequencies(&mut &src[..]);
        assert!(result.is_err_and(|e| e.kind() == io::ErrorKind::InvalidData));

        Ok(())
    }
}
