use std::io;

use super::{LOWER_BOUND, state_renormalize, state_step, write_alphabet, write_states};
use crate::{codecs::rans_nx16::ALPHABET_SIZE, io::writer::num::write_uint7};

const NORMALIZATION_BITS: u32 = 12;

type Frequencies = [u32; ALPHABET_SIZE];

pub(super) struct Context {
    frequencies: Frequencies,
}

pub(super) fn build_context(src: &[u8]) -> Context {
    let frequencies = build_frequencies(src);
    let normalized_frequencies = normalize_frequencies(&frequencies);

    Context {
        frequencies: normalized_frequencies,
    }
}

pub(super) fn write_context(dst: &mut Vec<u8>, ctx: &Context) -> io::Result<()> {
    let alphabet = build_alphabet(&ctx.frequencies);
    write_alphabet(dst, &alphabet)?;

    write_frequencies(dst, &ctx.frequencies)?;

    Ok(())
}

pub fn encode(src: &[u8], ctx: &Context, state_count: usize, dst: &mut Vec<u8>) -> io::Result<()> {
    let frequencies = &ctx.frequencies;
    let cumulative_frequencies = build_cumulative_frequencies(frequencies);

    let mut buf = Vec::new();
    let mut states = vec![LOWER_BOUND; state_count];

    for (i, &sym) in src.iter().enumerate().rev() {
        let state = &mut states[i % state_count];
        let j = usize::from(sym);
        let (f, g) = (frequencies[j], cumulative_frequencies[j]);
        *state = state_renormalize(*state, f, NORMALIZATION_BITS, &mut buf);
        *state = state_step(*state, f, g, NORMALIZATION_BITS);
    }

    write_states(dst, &states)?;
    dst.extend(buf.iter().rev());

    Ok(())
}

fn write_frequencies(dst: &mut Vec<u8>, frequencies: &Frequencies) -> io::Result<()> {
    for &f in frequencies {
        if f > 0 {
            write_uint7(dst, f)?;
        }
    }

    Ok(())
}

fn build_alphabet(frequencies: &Frequencies) -> [bool; ALPHABET_SIZE] {
    frequencies.map(|f| f > 0)
}

fn build_frequencies(src: &[u8]) -> Frequencies {
    let mut frequencies = [0; ALPHABET_SIZE];

    for &b in src {
        let i = usize::from(b);
        frequencies[i] += 1;
    }

    frequencies
}

pub(super) fn normalize_frequencies(frequencies: &Frequencies) -> Frequencies {
    // § 3.1 "Frequency tables" (2023-03-15): "...the final power of two used in the Order-0
    // (4096)..."
    const SCALING_FACTOR: u32 = 4096;

    let (max_index, sum) = describe_frequencies(frequencies);

    if sum == 0 {
        return [0; ALPHABET_SIZE];
    }

    let mut normalized_frequencies = [0; ALPHABET_SIZE];
    let mut normalized_sum = 0;

    for (&f, g) in frequencies.iter().zip(&mut normalized_frequencies) {
        if f == 0 {
            continue;
        }

        *g = (f * SCALING_FACTOR / sum).max(1);

        normalized_sum += *g;
    }

    if normalized_sum < SCALING_FACTOR {
        normalized_frequencies[max_index] += SCALING_FACTOR - normalized_sum;
    } else if normalized_sum > SCALING_FACTOR {
        normalized_frequencies[max_index] -= normalized_sum - SCALING_FACTOR;
    }

    normalized_frequencies
}

fn describe_frequencies(frequencies: &Frequencies) -> (usize, u32) {
    let mut max = u32::MIN;
    let mut max_index = 0;
    let mut sum = 0;

    for (i, &f) in frequencies.iter().enumerate() {
        if f >= max {
            max = f;
            max_index = i;
        }

        sum += f;
    }

    (max_index, sum)
}

pub(super) fn build_cumulative_frequencies(frequencies: &Frequencies) -> [u32; ALPHABET_SIZE] {
    let mut cumulative_frequencies = [0; ALPHABET_SIZE];
    let mut f = cumulative_frequencies[0];

    for (next_f, g) in cumulative_frequencies[1..].iter_mut().zip(frequencies) {
        *next_f = f + g;
        f = *next_f;
    }

    cumulative_frequencies
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_frequencies() -> io::Result<()> {
        let src = b"abracadabra";
        let frequencies = build_frequencies(src);

        let mut dst = Vec::new();
        write_frequencies(&mut dst, &frequencies)?;

        let expected = [
            0x05, // frequencies[b'a'] = 5
            0x02, // frequencies[b'b'] = 2
            0x01, // frequencies[b'c'] = 1
            0x01, // frequencies[b'd'] = 1
            0x02, // frequencies[b'r'] = 2
        ];

        assert_eq!(dst, expected);

        Ok(())
    }

    #[test]
    fn test_build_alphabet() {
        let src = b"abracadabra";
        let frequencies = build_frequencies(src);

        let mut expected = [false; ALPHABET_SIZE];

        for &sym in src {
            expected[usize::from(sym)] = true;
        }

        assert_eq!(build_alphabet(&frequencies), expected);
    }

    #[test]
    fn test_build_frequencies() {
        let src = b"abracadabra";

        let mut expected = [0; ALPHABET_SIZE];
        expected[usize::from(b'a')] = 5;
        expected[usize::from(b'b')] = 2;
        expected[usize::from(b'c')] = 1;
        expected[usize::from(b'd')] = 1;
        expected[usize::from(b'r')] = 2;

        assert_eq!(build_frequencies(src), expected);
    }

    #[test]
    fn test_normalize_frequencies() {
        let src = b"abracadabra";
        let frequencies = build_frequencies(src);

        let mut expected = [0; ALPHABET_SIZE];
        expected[usize::from(b'a')] = 1864;
        expected[usize::from(b'b')] = 744;
        expected[usize::from(b'c')] = 372;
        expected[usize::from(b'd')] = 372;
        expected[usize::from(b'r')] = 744;

        assert_eq!(normalize_frequencies(&frequencies), expected);
    }

    #[test]
    fn test_describe_frequencies() {
        let src = b"abracadabra";
        let frequencies = build_frequencies(src);
        let (max_index, sum) = describe_frequencies(&frequencies);
        assert_eq!(max_index, usize::from(b'a'));
        assert_eq!(sum, 11);
    }
}
