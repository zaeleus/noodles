use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::io::reader::num::read_uint7;

pub fn decode<R>(reader: &mut R, output: &mut [u8], n: u32) -> io::Result<()>
where
    R: Read,
{
    use super::{
        rans_advance_step_nx16, rans_get_cumulative_freq_nx16, rans_get_symbol_from_freq,
        rans_renorm_nx16,
    };

    let mut freqs = [0; 256];
    let mut cumulative_freqs = [0; 256];

    read_frequencies(reader, &mut freqs, &mut cumulative_freqs)?;

    let mut state = vec![0; n as usize];

    for s in &mut state {
        *s = reader.read_u32::<LittleEndian>()?;
    }

    for (i, b) in output.iter_mut().enumerate() {
        let j = i % (n as usize);

        let f = rans_get_cumulative_freq_nx16(state[j], 12);
        let s = rans_get_symbol_from_freq(&cumulative_freqs, f);

        *b = s;

        state[j] = rans_advance_step_nx16(
            state[j],
            cumulative_freqs[s as usize],
            freqs[s as usize],
            12,
        );

        state[j] = rans_renorm_nx16(reader, state[j])?;
    }

    Ok(())
}

pub fn normalize_frequencies(freqs: &mut [u32], bits: u32) {
    let mut total: u32 = freqs.iter().sum();

    if total == 0 || total == (1 << bits) {
        return;
    }

    let mut shift = 0;

    while total < (1 << bits) {
        total *= 2;
        shift += 1;
    }

    for freq in freqs {
        *freq <<= shift;
    }
}

fn read_frequencies<R>(
    reader: &mut R,
    freqs: &mut [u32],
    cumulative_freqs: &mut [u32],
) -> io::Result<()>
where
    R: Read,
{
    use super::read_alphabet;

    let alphabet = read_alphabet(reader)?;

    for i in 0..alphabet.len() {
        if alphabet[i] {
            freqs[i] = read_uint7(reader)?;
        }
    }

    normalize_frequencies(freqs, 12);

    cumulative_freqs[0] = 0;

    for i in 0..255 {
        cumulative_freqs[i + 1] = cumulative_freqs[i] + freqs[i];
    }

    Ok(())
}
