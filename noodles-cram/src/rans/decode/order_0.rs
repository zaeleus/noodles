use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::num::read_itf8;

use super::{rans_advance_step, rans_get_cumulative_freq, rans_get_symbol_from_freq, rans_renorm};

pub fn decode<R>(reader: &mut R, output: &mut [u8]) -> io::Result<()>
where
    R: Read,
{
    let mut freqs = vec![0; 256];
    let mut cumulative_freqs = vec![0; 256];

    read_frequencies_0(reader, &mut freqs, &mut cumulative_freqs)?;

    let mut state = [0; 4];
    reader.read_u32_into::<LittleEndian>(&mut state)?;

    let mut i = 0;

    while i < output.len() {
        for j in 0..4 {
            if i + j >= output.len() {
                return Ok(());
            }

            let f = rans_get_cumulative_freq(state[j]);
            let s = rans_get_symbol_from_freq(&cumulative_freqs, f);

            output[i + j] = s;

            state[j] = rans_advance_step(state[j], cumulative_freqs[s as usize], freqs[s as usize]);
            state[j] = rans_renorm(reader, state[j])?;
        }

        i += 4;
    }

    Ok(())
}

pub fn read_frequencies_0<R>(
    reader: &mut R,
    freqs: &mut [u32],
    cumulative_freqs: &mut [u32],
) -> io::Result<()>
where
    R: Read,
{
    let mut sym = reader.read_u8()?;
    let mut last_sym = sym;
    let mut rle = 0;

    loop {
        let f = read_itf8(reader)? as u32;

        freqs[sym as usize] = f;

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
