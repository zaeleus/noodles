use std::io;

use super::{
    TABLE_SIZE, build_table, decode_symbol, order_0, read_states, renormalize_state,
    renormalize_states,
};
use crate::{
    codecs::rans_4x8::{ALPHABET_SIZE, STATE_COUNT},
    io::reader::num::read_u8,
};

#[cfg(test)]
type Frequencies = [[u16; ALPHABET_SIZE]; ALPHABET_SIZE]; // F

/// The frequency rows read for each context. Most blocks have few contexts, so this avoids a
/// 256 x 256 table per block. `slots[context]` is the context's index in `rows`, if it has a row.
/// A context read again overwrites its row, as in a full table.
struct Rows {
    rows: Vec<[u16; ALPHABET_SIZE]>,
    slots: [Option<u8>; ALPHABET_SIZE],
}

impl Rows {
    fn set(&mut self, context: u8, row: [u16; ALPHABET_SIZE]) {
        let slot = &mut self.slots[usize::from(context)];

        match *slot {
            Some(i) => self.rows[usize::from(i)] = row,
            None => {
                // At most 256 rows, so this fits.
                *slot = Some(self.rows.len() as u8);
                self.rows.push(row);
            }
        }
    }
}

/// Decoding tables (see [`build_table`]) for every context with frequencies, stored
/// contiguously. Contexts without frequencies share an empty table (index 0), which decodes
/// like the all-zero row it represents. Quality scores can use around a hundred contexts, so the
/// compact per-context layout matters here.
struct Tables {
    symbols: Vec<[u8; TABLE_SIZE]>,
    frequencies: Vec<[u32; ALPHABET_SIZE]>,
    indices: [usize; ALPHABET_SIZE],
}

impl Tables {
    fn new(rows: &Rows) -> Self {
        let present = || {
            rows.slots.iter().enumerate().filter_map(|(i, &slot)| {
                let row = rows.rows.get(usize::from(slot?))?;
                row.iter().any(|&f| f != 0).then_some((i, row))
            })
        };

        let context_count = present().count();

        let mut tables = Self {
            symbols: vec![[0; TABLE_SIZE]; context_count + 1],
            frequencies: vec![[0; ALPHABET_SIZE]; context_count + 1],
            indices: [0; ALPHABET_SIZE],
        };

        build_table(
            &[0; ALPHABET_SIZE],
            &mut tables.symbols[0],
            &mut tables.frequencies[0],
        );

        for (index, (i, row)) in (1..).zip(present()) {
            build_table(
                row,
                &mut tables.symbols[index],
                &mut tables.frequencies[index],
            );
            tables.indices[i] = index;
        }

        tables
    }

    #[inline(always)]
    fn decode_symbol(&self, prev_sym: u8, state: &mut u32) -> u8 {
        let index = self.indices[usize::from(prev_sym)];
        decode_symbol(&self.symbols[index], &self.frequencies[index], state)
    }
}

pub fn decode(src: &mut &[u8], dst: &mut [u8]) -> io::Result<()> {
    let rows = read_frequency_rows(src)?;
    // Bounded totals keep every state step within u32, so the four steps can run ahead of the
    // ordered renormalization reads without overflowing on malformed tables.
    validate_frequencies(&rows)?;
    let tables = Tables::new(&rows);

    let mut states = read_states(src)?;
    let mut prev_syms = [0; STATE_COUNT];

    let [chunk_0, chunk_1, chunk_2, chunk_3, chunk_4] = split_chunks(dst);
    let chunks = chunk_0.iter_mut().zip(chunk_1).zip(chunk_2).zip(chunk_3);

    for (((d0, d1), d2), d3) in chunks {
        let dsts = [d0, d1, d2, d3];

        for (state, (prev_sym, d)) in states.iter_mut().zip(prev_syms.iter_mut().zip(dsts)) {
            let sym = tables.decode_symbol(*prev_sym, state);
            *d = sym;
            *prev_sym = sym;
        }

        // Keep the four independent state calculations ahead of the ordered byte reads. The
        // renormalization is branchless, since order-1 quality contexts make the per-state
        // renormalization branch unpredictable.
        renormalize_states(&mut states, src)?;
    }

    let mut state = states[3];
    let mut prev_sym = prev_syms[3];

    for d in chunk_4 {
        let sym = tables.decode_symbol(prev_sym, &mut state);
        *d = sym;
        state = renormalize_state(state, src)?;
        prev_sym = sym;
    }

    Ok(())
}

fn read_frequency_rows(src: &mut &[u8]) -> io::Result<Rows> {
    let mut rows = Rows {
        rows: Vec::new(),
        slots: [None; ALPHABET_SIZE],
    };

    let mut sym = read_u8(src)?;
    let mut prev_sym = sym;

    loop {
        let f = order_0::read_frequencies(src)?;
        rows.set(sym, f);

        sym = read_u8(src)?;

        if sym == 0 {
            break;
        }

        if sym - 1 == prev_sym {
            let len = read_u8(src)?;

            for _ in 0..len {
                let f = order_0::read_frequencies(src)?;
                rows.set(sym, f);
                // A malformed run can pass 255; wrap (as release builds did) instead of panicking.
                sym = sym.wrapping_add(1);
            }
        }

        prev_sym = sym;
    }

    Ok(rows)
}

#[cfg(test)]
fn read_frequencies(src: &mut &[u8]) -> io::Result<Frequencies> {
    let rows = read_frequency_rows(src)?;
    let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    for (row, &slot) in frequencies.iter_mut().zip(&rows.slots) {
        if let Some(i) = slot {
            *row = rows.rows[usize::from(i)];
        }
    }

    Ok(frequencies)
}

fn validate_frequencies(rows: &Rows) -> io::Result<()> {
    const TOTAL_FREQUENCY: u32 = 4096;

    // Only each context's last row counts, as in a full table.
    for row in &rows.rows {
        let total: u32 = row.iter().copied().map(u32::from).sum();

        if total > TOTAL_FREQUENCY {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid rANS order-1 frequency table",
            ));
        }
    }

    Ok(())
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
    fn test_decode_rejects_oversized_frequency_tables() {
        // Context '\0': 'a' = 4095, 'c' = 2 (total 4097), with no states following.
        let src = [0x00, b'a', 0x8f, 0xff, b'c', 0x02, 0x00, 0x00];
        let mut dst = [0; 8];
        assert_eq!(
            decode(&mut &src[..], &mut dst).unwrap_err().kind(),
            io::ErrorKind::InvalidData
        );
    }

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
