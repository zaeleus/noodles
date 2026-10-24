mod header;
mod order_0;
mod order_1;

use std::io;

use self::header::read_header;
use super::{ALPHABET_SIZE, LOWER_BOUND, Order, STATE_COUNT};
use crate::io::reader::num::{read_u8, read_u32_le};

pub fn decode(mut src: &[u8]) -> io::Result<Vec<u8>> {
    let (order, _, uncompressed_size) = read_header(&mut src)?;

    let mut dst = vec![0; uncompressed_size];

    match order {
        Order::Zero => order_0::decode(&mut src, &mut dst)?,
        Order::One => order_1::decode(&mut src, &mut dst)?,
    }

    Ok(dst)
}

fn read_states(src: &mut &[u8]) -> io::Result<[u32; STATE_COUNT]> {
    let mut states = [0; STATE_COUNT];

    for state in &mut states {
        *state = read_u32_le(src)?;
    }

    Ok(states)
}

/// The rANS 4x8 frequency precision (TF_SHIFT).
const NORMALIZATION_BITS: u32 = 12;

/// The number of cumulative frequency slots in a table.
const TABLE_SIZE: usize = 1 << NORMALIZATION_BITS;

/// Builds a compact decoding table from one frequency row whose total is at most 4096.
///
/// `symbols` maps each cumulative frequency slot to its symbol, and `frequencies` is indexed by
/// symbol and packs `frequency << 16 | cumulative_frequency`. (rANS Nx16 instead packs per slot;
/// see its `SymbolTable`.) At 5 KiB per table, this stays cache-resident and cheap to build for
/// small blocks.
fn build_table(
    row: &[u16; ALPHABET_SIZE],
    symbols: &mut [u8; TABLE_SIZE],
    frequencies: &mut [u32; ALPHABET_SIZE],
) {
    let mut cumulative_frequency = 0;

    for ((sym, &f), packed) in (0..=u8::MAX).zip(row).zip(frequencies.iter_mut()) {
        let next_cumulative_frequency = cumulative_frequency + usize::from(f);
        symbols[cumulative_frequency..next_cumulative_frequency].fill(sym);
        *packed = (u32::from(f) << 16) | cumulative_frequency as u32;
        cumulative_frequency = next_cumulative_frequency;
    }

    // Slots at or past the total frequency resolve to the last symbol, as with a search for the
    // last cumulative frequency not greater than the slot.
    symbols[cumulative_frequency..].fill(u8::MAX);
}

/// Decodes a symbol with a table from [`build_table`] and advances the state, before
/// renormalization.
#[inline(always)]
fn decode_symbol(
    symbols: &[u8; TABLE_SIZE],
    frequencies: &[u32; ALPHABET_SIZE],
    state: &mut u32,
) -> u8 {
    let f = (*state & 0x0fff) as usize;
    let sym = symbols[f];
    let packed = frequencies[usize::from(sym)];
    let (frequency, cumulative_frequency) = (packed >> 16, packed & 0xffff);

    // Valid tables cannot overflow here; wrapping keeps malformed initial states panic-free.
    *state = frequency
        .wrapping_mul(*state >> NORMALIZATION_BITS)
        .wrapping_add(f as u32)
        .wrapping_sub(cumulative_frequency);

    sym
}

/// Renormalizes one state, reading bytes from `src` as needed.
///
/// A step with a nonzero symbol frequency leaves the state at or above 2^11, so at most two
/// bytes are needed. Those are applied without branching on whether the state needs them.
/// A larger shortfall (only possible with malformed input) continues byte by byte, so bytes are
/// consumed in the same order, and errors are the same, as a plain byte loop.
///
/// Each step sets `r` to 1 if the state needs a byte and 0 otherwise. Shifting by `8 * r` and
/// masking the byte with `r.wrapping_neg()` (all ones or all zeros) then either applies the step
/// or leaves the state unchanged.
#[inline(always)]
fn renormalize_state(mut s: u32, src: &mut &[u8]) -> io::Result<u32> {
    if let [b0, b1, ..] = **src {
        let r0 = u32::from(s < LOWER_BOUND);
        s = (s << (8 * r0)) | (u32::from(b0) & r0.wrapping_neg());

        let r1 = u32::from(s < LOWER_BOUND);
        s = (s << (8 * r1)) | (u32::from(b1) & r1.wrapping_neg());

        *src = &src[(r0 + r1) as usize..];

        if s >= LOWER_BOUND {
            return Ok(s);
        }
    }

    state_renormalize(s, src)
}

/// Renormalizes the states in order, after all of their independent steps.
#[inline(always)]
fn renormalize_states(states: &mut [u32; STATE_COUNT], src: &mut &[u8]) -> io::Result<()> {
    // Each state takes at most two bytes from a valid stream, so one length check covers all
    // four. The steps are as in `renormalize_state`. State `i` starts at most at byte `2 * i`,
    // so `k` is below 8 at every load and `k & 7 == k`; the mask only elides the bounds checks.
    if let Some(window) = src.first_chunk::<{ 2 * STATE_COUNT }>() {
        let mut k = 0;

        for (i, state) in states.iter_mut().enumerate() {
            let mut s = *state;

            let r0 = u32::from(s < LOWER_BOUND);
            s = (s << (8 * r0)) | (u32::from(window[k & 7]) & r0.wrapping_neg());
            k += r0 as usize;

            // A second byte is rarely needed, so this branch predicts well.
            if s < LOWER_BOUND {
                s = (s << 8) | u32::from(window[k & 7]);
                k += 1;
            }

            *state = s;

            if s < LOWER_BOUND {
                // Only malformed streams need a third byte. Finish this state and the rest
                // byte by byte, in order.
                *src = &src[k..];
                *state = state_renormalize(s, src)?;

                for state in &mut states[i + 1..] {
                    *state = renormalize_state(*state, src)?;
                }

                return Ok(());
            }
        }

        *src = &src[k..];
        return Ok(());
    }

    for state in states {
        *state = renormalize_state(*state, src)?;
    }

    Ok(())
}

fn state_renormalize(mut s: u32, src: &mut &[u8]) -> io::Result<u32> {
    while s < LOWER_BOUND {
        let b = read_u8(src).map(u32::from)?;
        s = (s << 8) | b;
    }

    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_states() -> io::Result<()> {
        let src = [
            0x05, 0x00, 0x00, 0x00, // states[0] = 5
            0x08, 0x00, 0x00, 0x00, // states[1] = 8
            0x0d, 0x00, 0x00, 0x00, // states[2] = 13
            0x15, 0x00, 0x00, 0x00, // states[3] = 21
        ];

        assert_eq!(read_states(&mut &src[..])?, [5, 8, 13, 21]);

        Ok(())
    }

    #[test]
    fn test_renormalize_state() {
        // Up to two bytes in the fast path, then byte by byte (only for malformed streams).
        for (state, src, expected) in [
            (LOWER_BOUND, &[1][..], Some((LOWER_BOUND, &[1][..]))),
            (0x1234, &[0xab, 0xcd, 0xef], Some((0x1234abcd, &[0xef][..]))),
            (1, &[2, 3, 4, 5], Some((0x01020304, &[5][..]))),
            (1, &[2], None),
        ] {
            let mut rest = src;
            let actual = renormalize_state(state, &mut rest).ok().map(|s| (s, rest));
            assert_eq!(actual, expected, "state={state:#x}");
        }
    }

    #[test]
    fn test_decode_with_order_0() -> io::Result<()> {
        let expected = b"noodles";

        let src = [
            0x00, 0x25, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x64, 0x82, 0x49, 0x65, 0x00,
            0x82, 0x49, 0x6c, 0x82, 0x49, 0x6e, 0x82, 0x49, 0x6f, 0x00, 0x84, 0x92, 0x73, 0x82,
            0x49, 0x00, 0xe2, 0x06, 0x83, 0x18, 0x74, 0x7b, 0x41, 0x0c, 0x2b, 0xa9, 0x41, 0x0c,
            0x25, 0x31, 0x80, 0x03,
        ];

        let actual = decode(&src)?;
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_decode_with_order_1() -> io::Result<()> {
        let expected = b"noodles";

        let src = [
            0x01, 0x3b, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x00, 0x64, 0x84, 0x00, 0x6e,
            0x84, 0x00, 0x6f, 0x00, 0x87, 0xff, 0x00, 0x64, 0x6c, 0x8f, 0xff, 0x00, 0x65, 0x00,
            0x73, 0x8f, 0xff, 0x00, 0x6c, 0x65, 0x8f, 0xff, 0x00, 0x6e, 0x6f, 0x8f, 0xff, 0x00,
            0x6f, 0x00, 0x64, 0x87, 0xff, 0x6f, 0x88, 0x00, 0x00, 0x00, 0x00, 0x04, 0x00, 0x02,
            0x02, 0x28, 0x00, 0x01, 0x02, 0x28, 0x00, 0x01, 0x02, 0x60, 0x00, 0x02,
        ];

        let actual = decode(&src)?;
        assert_eq!(actual, expected);

        Ok(())
    }
}
