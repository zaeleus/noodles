use std::io::{self, Read};

use crate::io::reader::num::{read_u8, read_u32_le, read_uint7, read_uint7_as};

pub fn decode<R>(reader: &mut R, output: &mut [u8], state_count: usize) -> io::Result<()>
where
    R: Read,
{
    use super::{
        rans_advance_step_nx16, rans_get_cumulative_freq_nx16, rans_get_symbol_from_freq,
        rans_renorm_nx16,
    };

    let mut freqs = vec![vec![0; 256]; 256];
    let mut cumulative_freqs = vec![vec![0; 256]; 256];

    let bits = read_frequencies(reader, &mut freqs, &mut cumulative_freqs)?;

    let mut state = vec![0; state_count];

    for s in &mut state {
        *s = read_u32_le(reader)?;
    }

    let mut i = 0;
    let mut last_syms = vec![0; state.len()];

    while i < output.len() / state_count {
        for j in 0..state_count {
            let f = rans_get_cumulative_freq_nx16(state[j], bits);
            let s = rans_get_symbol_from_freq(&cumulative_freqs[last_syms[j]], f);

            output[i + j * (output.len() / state_count)] = s;

            state[j] = rans_advance_step_nx16(
                state[j],
                cumulative_freqs[last_syms[j]][usize::from(s)],
                freqs[last_syms[j]][usize::from(s)],
                bits,
            );

            state[j] = rans_renorm_nx16(reader, state[j])?;

            last_syms[j] = usize::from(s);
        }

        i += 1;
    }

    i *= state_count;
    let m = state_count - 1;

    while i < output.len() {
        let f = rans_get_cumulative_freq_nx16(state[m], bits);
        let s = rans_get_symbol_from_freq(&cumulative_freqs[last_syms[m]], f);

        output[i] = s;

        state[m] = rans_advance_step_nx16(
            state[m],
            cumulative_freqs[last_syms[m]][usize::from(s)],
            freqs[last_syms[m]][usize::from(s)],
            bits,
        );

        state[m] = rans_renorm_nx16(reader, state[m])?;

        last_syms[m] = usize::from(s);

        i += 1;
    }

    Ok(())
}

fn read_frequencies<R>(
    reader: &mut R,
    freqs: &mut [Vec<u32>],
    cumulative_freqs: &mut [Vec<u32>],
) -> io::Result<u32>
where
    R: Read,
{
    use super::order_0;

    let comp = read_u8(reader)?;
    let bits = u32::from(comp >> 4);

    if comp & 0x01 != 0 {
        let u_size = read_uint7_as(reader)?;
        let c_size = read_uint7_as(reader)?;

        let mut c_data = vec![0; c_size];
        reader.read_exact(&mut c_data)?;

        let mut c_data_reader = &c_data[..];
        let mut u_data = vec![0; u_size];
        order_0::decode(&mut c_data_reader, &mut u_data, 4)?;

        let mut u_data_reader = &u_data[..];
        read_frequencies_inner(&mut u_data_reader, freqs, cumulative_freqs, bits)?;
    } else {
        read_frequencies_inner(reader, freqs, cumulative_freqs, bits)?;
    }

    Ok(bits)
}

fn read_frequencies_inner<R>(
    reader: &mut R,
    freqs: &mut [Vec<u32>],
    cumulative_freqs: &mut [Vec<u32>],
    bits: u32,
) -> io::Result<()>
where
    R: Read,
{
    use super::{order_0, read_alphabet};

    let alphabet = read_alphabet(reader)?;

    for (i, a) in alphabet.iter().enumerate() {
        if !a {
            continue;
        }

        let mut run = 0;

        for (j, b) in alphabet.iter().enumerate() {
            if !b {
                continue;
            }

            if run > 0 {
                run -= 1;
            } else {
                let f = read_uint7(reader)?;

                freqs[i][j] = f;

                if f == 0 {
                    run = read_u8(reader)?;
                }
            }
        }

        order_0::normalize_frequencies(&mut freqs[i], bits);

        cumulative_freqs[i][0] = 0;

        for j in 0..255 {
            cumulative_freqs[i][j + 1] = cumulative_freqs[i][j] + freqs[i][j];
        }
    }

    Ok(())
}
