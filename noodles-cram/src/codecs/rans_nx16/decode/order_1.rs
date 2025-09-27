use std::io::{self, Read};

use super::read_states;
use crate::io::reader::num::{read_u8, read_uint7, read_uint7_as};

pub fn decode(src: &mut &[u8], dst: &mut [u8], state_count: usize) -> io::Result<()> {
    use super::{
        rans_advance_step_nx16, rans_get_cumulative_freq_nx16, rans_get_symbol_from_freq,
        rans_renorm_nx16,
    };

    let mut freqs = vec![vec![0; 256]; 256];
    let mut cumulative_freqs = vec![vec![0; 256]; 256];

    let bits = read_frequencies(src, &mut freqs, &mut cumulative_freqs)?;

    let mut states = read_states(src, state_count)?;

    let mut i = 0;
    let mut last_syms = vec![0; states.len()];

    while i < dst.len() / state_count {
        for j in 0..state_count {
            let f = rans_get_cumulative_freq_nx16(states[j], bits);
            let s = rans_get_symbol_from_freq(&cumulative_freqs[last_syms[j]], f);

            dst[i + j * (dst.len() / state_count)] = s;

            states[j] = rans_advance_step_nx16(
                states[j],
                cumulative_freqs[last_syms[j]][usize::from(s)],
                freqs[last_syms[j]][usize::from(s)],
                bits,
            );

            states[j] = rans_renorm_nx16(src, states[j])?;

            last_syms[j] = usize::from(s);
        }

        i += 1;
    }

    i *= state_count;
    let m = state_count - 1;

    while i < dst.len() {
        let f = rans_get_cumulative_freq_nx16(states[m], bits);
        let s = rans_get_symbol_from_freq(&cumulative_freqs[last_syms[m]], f);

        dst[i] = s;

        states[m] = rans_advance_step_nx16(
            states[m],
            cumulative_freqs[last_syms[m]][usize::from(s)],
            freqs[last_syms[m]][usize::from(s)],
            bits,
        );

        states[m] = rans_renorm_nx16(src, states[m])?;

        last_syms[m] = usize::from(s);

        i += 1;
    }

    Ok(())
}

fn read_frequencies(
    src: &mut &[u8],
    freqs: &mut [Vec<u32>],
    cumulative_freqs: &mut [Vec<u32>],
) -> io::Result<u32> {
    use super::order_0;

    let comp = read_u8(src)?;
    let bits = u32::from(comp >> 4);

    if comp & 0x01 != 0 {
        let u_size = read_uint7_as(src)?;
        let c_size = read_uint7_as(src)?;

        let mut c_data = vec![0; c_size];
        src.read_exact(&mut c_data)?;

        let mut c_data_reader = &c_data[..];
        let mut u_data = vec![0; u_size];
        order_0::decode(&mut c_data_reader, &mut u_data, 4)?;

        let mut u_data_reader = &u_data[..];
        read_frequencies_inner(&mut u_data_reader, freqs, cumulative_freqs, bits)?;
    } else {
        read_frequencies_inner(src, freqs, cumulative_freqs, bits)?;
    }

    Ok(bits)
}

fn read_frequencies_inner(
    src: &mut &[u8],
    freqs: &mut [Vec<u32>],
    cumulative_freqs: &mut [Vec<u32>],
    bits: u32,
) -> io::Result<()> {
    use super::{order_0, read_alphabet};

    let alphabet = read_alphabet(src)?;

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
                let f = read_uint7(src)?;

                freqs[i][j] = f;

                if f == 0 {
                    run = read_u8(src)?;
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
