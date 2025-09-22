use std::io;

use crate::io::reader::num::{read_u32_le, read_uint7};

pub fn decode(src: &mut &[u8], dst: &mut [u8], state_count: usize) -> io::Result<()> {
    use super::{
        rans_advance_step_nx16, rans_get_cumulative_freq_nx16, rans_get_symbol_from_freq,
        rans_renorm_nx16,
    };

    let mut freqs = [0; 256];
    let mut cumulative_freqs = [0; 256];

    read_frequencies(src, &mut freqs, &mut cumulative_freqs)?;

    let mut state = vec![0; state_count];

    for s in &mut state {
        *s = read_u32_le(src)?;
    }

    for (i, d) in dst.iter_mut().enumerate() {
        let j = i % state_count;

        let f = rans_get_cumulative_freq_nx16(state[j], 12);
        let s = rans_get_symbol_from_freq(&cumulative_freqs, f);

        *d = s;

        state[j] = rans_advance_step_nx16(
            state[j],
            cumulative_freqs[s as usize],
            freqs[s as usize],
            12,
        );

        state[j] = rans_renorm_nx16(src, state[j])?;
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

fn read_frequencies(
    src: &mut &[u8],
    freqs: &mut [u32],
    cumulative_freqs: &mut [u32],
) -> io::Result<()> {
    use super::read_alphabet;

    let alphabet = read_alphabet(src)?;

    for i in 0..alphabet.len() {
        if alphabet[i] {
            freqs[i] = read_uint7(src)?;
        }
    }

    normalize_frequencies(freqs, 12);

    cumulative_freqs[0] = 0;

    for i in 0..255 {
        cumulative_freqs[i + 1] = cumulative_freqs[i] + freqs[i];
    }

    Ok(())
}
