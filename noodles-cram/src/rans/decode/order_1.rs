use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use super::{
    order_0, rans_advance_step, rans_get_cumulative_freq, rans_get_symbol_from_freq, rans_renorm,
};

pub fn decode<R>(reader: &mut R, output: &mut [u8]) -> io::Result<()>
where
    R: Read,
{
    let mut freqs = vec![vec![0; 256]; 256];
    let mut cumulative_freqs = vec![vec![0; 256]; 256];

    read_frequencies_1(reader, &mut freqs, &mut cumulative_freqs)?;

    let mut state = [0; 4];
    reader.read_u32_into::<LittleEndian>(&mut state)?;

    let mut i = 0;
    let mut last_syms = [0; 4];

    while i < output.len() / 4 {
        for j in 0..4 {
            let f = rans_get_cumulative_freq(state[j]);
            let s = rans_get_symbol_from_freq(&cumulative_freqs[last_syms[j] as usize], f);

            output[i + j * (output.len() / 4)] = s;

            state[j] = rans_advance_step(
                state[j],
                cumulative_freqs[last_syms[j] as usize][s as usize],
                freqs[last_syms[j] as usize][s as usize],
            );
            state[j] = rans_renorm(reader, state[j])?;

            last_syms[j] = s;
        }

        i += 1;
    }

    i *= 4;

    while i < output.len() {
        let f = rans_get_cumulative_freq(state[3]);
        let s = rans_get_symbol_from_freq(&cumulative_freqs[last_syms[3] as usize], f);

        output[i] = s;

        state[3] = rans_advance_step(
            state[3],
            cumulative_freqs[last_syms[3] as usize][s as usize],
            freqs[last_syms[3] as usize][s as usize],
        );
        state[3] = rans_renorm(reader, state[3])?;

        last_syms[3] = s;

        i += 1;
    }

    Ok(())
}

fn read_frequencies_1<R>(
    reader: &mut R,
    freqs: &mut Vec<Vec<u32>>,
    cumulative_freqs: &mut Vec<Vec<u32>>,
) -> io::Result<()>
where
    R: Read,
{
    let mut sym = reader.read_u8()?;
    let mut last_sym = sym;
    let mut rle = 0;

    loop {
        order_0::read_frequencies_0(
            reader,
            &mut freqs[sym as usize],
            &mut cumulative_freqs[sym as usize],
        )?;

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

    Ok(())
}
