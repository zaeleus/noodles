use std::io;

use super::{
    NORMALIZATION_BITS, TABLE_SIZE, build_table, decode_symbol, read_states, state_renormalize,
};
use crate::{
    codecs::rans_4x8::{ALPHABET_SIZE, STATE_COUNT},
    io::reader::num::{read_itf8_as, read_u8},
};

type Frequencies = [u16; ALPHABET_SIZE]; // F

pub fn decode(src: &mut &[u8], dst: &mut [u8]) -> io::Result<()> {
    let frequencies = read_frequencies(src)?;
    validate_frequencies(&frequencies)?;

    let mut symbols = [0; TABLE_SIZE];
    let mut packed_frequencies = [0; ALPHABET_SIZE];
    build_table(&frequencies, &mut symbols, &mut packed_frequencies);

    let mut states = read_states(src)?;
    let (chunks, remainder) = dst.as_chunks_mut::<STATE_COUNT>();

    for chunk in chunks {
        for (d, state) in chunk.iter_mut().zip(&mut states) {
            *d = decode_symbol(&symbols, &packed_frequencies, state);
        }

        // The state steps are independent, so they run ahead of the ordered renormalization reads.
        // Unlike order-1, the plain byte loop is faster here than the branchless
        // renormalization.
        for state in &mut states {
            *state = state_renormalize(*state, src)?;
        }
    }

    for (d, state) in remainder.iter_mut().zip(&mut states) {
        *d = decode_symbol(&symbols, &packed_frequencies, state);
        *state = state_renormalize(*state, src)?;
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
                // A malformed run can pass 255; wrap (as release builds did) instead of panicking.
                sym = sym.wrapping_add(1);
            }
        }

        prev_sym = sym;
    }

    Ok(frequencies)
}

/// Rejects tables whose total exceeds 4096, which cannot come from a valid encoder and would
/// overflow the cumulative frequencies.
fn validate_frequencies(frequencies: &Frequencies) -> io::Result<()> {
    const TOTAL_FREQUENCY: u32 = 1 << NORMALIZATION_BITS;

    let total: u32 = frequencies.iter().copied().map(u32::from).sum();

    if total > TOTAL_FREQUENCY {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid rANS order-0 frequency table",
        ))
    } else {
        Ok(())
    }
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
}
