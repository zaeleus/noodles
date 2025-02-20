use std::io::{self, Read};

use byteorder::ReadBytesExt;

use crate::io::reader::num::read_itf8;

use super::{rans_advance_step, rans_get_cumulative_freq, rans_renorm, read_states};

pub fn decode<R>(reader: &mut R, dst: &mut [u8]) -> io::Result<()>
where
    R: Read,
{
    let mut freqs = [0; 256];
    let mut cumulative_freqs = [0; 256];

    read_frequencies_0(reader, &mut freqs, &mut cumulative_freqs)?;

    let cumulative_freqs_symbols_table = build_cumulative_freqs_symbols_table_0(&cumulative_freqs);

    let mut states = read_states(reader)?;

    for chunk in dst.chunks_mut(states.len()) {
        for (d, state) in chunk.iter_mut().zip(states.iter_mut()) {
            let f = rans_get_cumulative_freq(*state);
            let s = cumulative_freqs_symbols_table[f as usize];

            *d = s;

            *state = rans_advance_step(
                *state,
                cumulative_freqs[usize::from(s)],
                freqs[usize::from(s)],
            );

            *state = rans_renorm(reader, *state)?;
        }
    }

    Ok(())
}

pub fn read_frequencies_0<R>(
    reader: &mut R,
    freqs: &mut [u32; 256],
    cumulative_freqs: &mut [u32; 256],
) -> io::Result<()>
where
    R: Read,
{
    let mut sym = reader.read_u8()?;
    let mut last_sym = sym;
    let mut rle = 0;

    loop {
        let f = read_itf8(reader)? as u32;

        freqs[usize::from(sym)] = f;

        if rle > 0 {
            rle -= 1;
            sym += 1;
        } else {
            sym = reader.read_u8()?;

            if last_sym < 255 && sym == last_sym + 1 {
                rle = reader.read_u8()?;
            }
        }

        last_sym = sym;

        if sym == 0 {
            break;
        }
    }

    cumulative_freqs[0] = 0;

    for i in 0..255 {
        cumulative_freqs[i + 1] = cumulative_freqs[i] + freqs[i];
    }

    Ok(())
}

pub fn build_cumulative_freqs_symbols_table_0(cumulative_freqs: &[u32; 256]) -> [u8; 4096] {
    let mut table = [0; 4096];
    let mut sym = 0;

    for (freq, cumulative_freq) in table.iter_mut().enumerate() {
        let freq = freq as u32;

        while sym < 255 && freq >= cumulative_freqs[usize::from(sym + 1)] {
            sym += 1;
        }

        *cumulative_freq = sym;
    }

    table
}
