pub mod bit_pack;
mod order_0;
mod order_1;
mod rle;
mod stripe;

use std::io;

use super::{ALPHABET_SIZE, Flags};
use crate::io::reader::num::{read_u8, read_u32_le, read_uint7_as};

pub fn decode(mut src: &[u8], mut uncompressed_size: usize) -> io::Result<Vec<u8>> {
    let flags = read_flags(&mut src)?;

    let state_count = flags.state_count();

    if flags.has_uncompressed_size() {
        uncompressed_size = read_uncompressed_size(&mut src)?;
    }

    if flags.is_striped() {
        return stripe::decode(&mut src, uncompressed_size);
    }

    let bit_pack_context = if flags.is_bit_packed() {
        let (ctx, len) = bit_pack::read_context(&mut src, uncompressed_size)?;
        uncompressed_size = len;
        Some(ctx)
    } else {
        None
    };

    let rle_context = if flags.is_rle() {
        let (ctx, len) = rle::read_context(&mut src, state_count, uncompressed_size)?;
        uncompressed_size = len;
        Some(ctx)
    } else {
        None
    };

    let mut dst = vec![0; uncompressed_size];

    if flags.is_uncompressed() {
        dst.copy_from_slice(src);
    } else if flags.order() == 0 {
        order_0::decode(&mut src, &mut dst, state_count)?;
    } else {
        order_1::decode(&mut src, &mut dst, state_count)?;
    };

    if let Some(ctx) = rle_context {
        dst = rle::decode(&dst, &ctx)?;
    }

    if let Some(ctx) = bit_pack_context {
        dst = bit_pack::decode(&dst, &ctx)?;
    }

    Ok(dst)
}

fn read_flags(src: &mut &[u8]) -> io::Result<Flags> {
    read_u8(src).map(Flags::from)
}

fn read_uncompressed_size(src: &mut &[u8]) -> io::Result<usize> {
    read_uint7_as(src)
}

fn read_alphabet(src: &mut &[u8]) -> io::Result<[bool; ALPHABET_SIZE]> {
    const NUL: u8 = 0x00;

    let mut alphabet = [false; ALPHABET_SIZE];

    let mut sym = read_u8(src)?;
    let mut prev_sym = sym;

    loop {
        alphabet[usize::from(sym)] = true;

        sym = read_u8(src)?;

        if sym == NUL {
            break;
        }

        if sym - 1 == prev_sym {
            let len = read_u8(src)?;

            for _ in 0..len {
                alphabet[usize::from(sym)] = true;
                // A malformed run can pass 255; wrap (as release builds did) instead of panicking.
                sym = sym.wrapping_add(1);
            }
        }

        prev_sym = sym;
    }

    Ok(alphabet)
}

fn read_states<const N: usize>(src: &mut &[u8]) -> io::Result<[u32; N]> {
    let mut states = [0; N];

    for state in &mut states {
        *state = read_u32_le(src)?;
    }

    Ok(states)
}

/// A decoding table for one frequency distribution.
///
/// Both arrays are indexed by cumulative frequency slot `f` (`1 << bits` entries). `symbols[f]` is
/// the symbol whose cumulative frequency range contains `f`, and `steps[f]` packs
/// `frequency << 16 | (f - cumulative_frequency)` for that symbol. (rANS 4x8 instead packs
/// `frequency << 16 | cumulative_frequency` per symbol.) This replaces a per-symbol search over
/// the cumulative frequencies, and the two lookups are independent of each other.
struct SymbolTable {
    symbols: Box<[u8]>,
    steps: Box<[u32]>,
}

impl SymbolTable {
    fn new(frequencies: &[u32; ALPHABET_SIZE], bits: u32) -> io::Result<Self> {
        let size = 1 << bits;
        let mut symbols = vec![0; size].into_boxed_slice();
        let mut steps = vec![0; size].into_boxed_slice();
        build_symbol_table(frequencies, bits, &mut symbols, &mut steps)?;
        Ok(Self { symbols, steps })
    }

    /// Decodes a symbol and advances the state, before renormalization.
    #[inline(always)]
    fn decode_symbol(&self, state: &mut u32, bits: u32) -> u8 {
        let f = (*state & ((1 << bits) - 1)) as usize;
        let sym = self.symbols[f];
        *state = apply_step(*state, self.steps[f], bits);
        sym
    }
}

/// Fills `symbols` and `steps` (each `1 << bits` long) as described for [`SymbolTable`].
fn build_symbol_table(
    frequencies: &[u32; ALPHABET_SIZE],
    bits: u32,
    symbols: &mut [u8],
    steps: &mut [u32],
) -> io::Result<()> {
    let total_frequency = 1 << bits;
    let mut cumulative_frequency: u32 = 0;

    for (sym, &f) in (0..=u8::MAX).zip(frequencies) {
        let next_cumulative_frequency = cumulative_frequency
            .checked_add(f)
            .filter(|&n| n <= total_frequency)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid rANS Nx16 frequency table",
                )
            })?;

        let range = cumulative_frequency as usize..next_cumulative_frequency as usize;
        symbols[range.clone()].fill(sym);

        // `f` and `offset` are at most `1 << bits`, and `bits` is at most 12 (order 0 uses 12;
        // order 1 rejects more), so both fit in 16 bits.
        for (offset, step) in (0..).zip(&mut steps[range]) {
            *step = (f << 16) | offset;
        }

        cumulative_frequency = next_cumulative_frequency;
    }

    // Slots at or past the total frequency resolve to the last symbol, as with a search for the
    // last cumulative frequency not greater than the slot.
    let last_frequency = frequencies[usize::from(u8::MAX)];
    let last_cumulative_frequency = cumulative_frequency - last_frequency;

    for (slot, (sym, step)) in (cumulative_frequency..).zip(
        symbols
            .iter_mut()
            .zip(steps.iter_mut())
            .skip(cumulative_frequency as usize),
    ) {
        *sym = u8::MAX;
        *step = (last_frequency << 16) | (slot - last_cumulative_frequency);
    }

    Ok(())
}

/// Advances a state with a packed step from a [`SymbolTable`]-style table.
#[inline(always)]
fn apply_step(state: u32, step: u32, bits: u32) -> u32 {
    // Valid tables cannot overflow here; wrapping keeps malformed initial states panic-free.
    (step >> 16)
        .wrapping_mul(state >> bits)
        .wrapping_add(step & 0xffff)
}

/// Renormalizes every state in order, reading from `src` as needed.
#[inline(always)]
fn renormalize_states<const N: usize>(states: &mut [u32; N], src: &mut &[u8]) -> io::Result<()> {
    if let Some((window, _)) = src.split_at_checked(2 * N) {
        // Fast path: enough input for every state to read, so no per-state EOF handling. Whether
        // a state renormalizes is data dependent and poorly predicted, so this is branchless:
        // every state loads the next two bytes, and only the states that need them use them.
        // `n <= 2 * i` before state `i`, so the loads stay within `window`.
        let mut n = 0;

        for state in states.iter_mut() {
            let renormalize = u32::from(*state < (1 << 15));
            let b = u32::from(u16::from_le_bytes([window[n], window[n + 1]]));
            *state = (*state << (16 * renormalize)) | (b & renormalize.wrapping_neg());
            n += 2 * renormalize as usize;
        }

        *src = &src[n..];
    } else {
        for state in states.iter_mut() {
            *state = state_renormalize(*state, src)?;
        }
    }

    Ok(())
}

fn state_renormalize(mut s: u32, src: &mut &[u8]) -> io::Result<u32> {
    if s < (1 << 15) {
        s = (s << 16) + read_u16_le(src).map(u32::from)?;
    }

    Ok(s)
}

fn read_u16_le(src: &mut &[u8]) -> io::Result<u16> {
    let (buf, rest) = src
        .split_first_chunk()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Ok(u16::from_le_bytes(*buf))
}

fn split_off<'a>(src: &mut &'a [u8], len: usize) -> io::Result<&'a [u8]> {
    let (buf, rest) = src
        .split_at_checked(len)
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_order_0() -> io::Result<()> {
        let src = [
            0x00, // flags = {empty}
            0x07, // uncompressed len = 7
            0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x01, 0x01, 0x01, 0x01, 0x03,
            0x01, 0x00, 0x26, 0x20, 0x00, 0x00, 0xb8, 0x0a, 0x00, 0x00, 0xd8, 0x0a, 0x00, 0x00,
            0x00, 0x04, 0x00,
        ];

        assert_eq!(decode(&src, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_decode_order_1() -> io::Result<()> {
        let src = [
            0x01, // flags = ORDER
            0x4d, // uncompressed len = 77
            0xa0, 0x00, 0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x00, 0x00, 0x01,
            0x01, 0x00, 0x00, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x01, 0x00,
            0x02, 0x00, 0x01, 0x0f, 0x00, 0x02, 0x01, 0x00, 0x01, 0x01, 0x0f, 0x00, 0x02, 0x00,
            0x03, 0x0f, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x02, 0x0f, 0x00, 0x00, 0x00,
            0x05, 0x10, 0x80, 0x72, 0x60, 0x00, 0x80, 0x8b, 0x5f, 0x00, 0xc0, 0xb0, 0x60, 0x00,
            0x40, 0x49, 0x39, 0x00,
        ];

        assert_eq!(
            decode(&src, 0)?,
            b"nnnnnnnnnnnnooooooooooooooooddddddddddddddllllllllllllllleeeeeeeeeessssssssss"
        );

        Ok(())
    }

    #[test]
    fn test_decode_stripe() -> io::Result<()> {
        let src = [
            0x08, // flags = STRIPE
            0x07, // uncompressed len = 7
            0x04, 0x17, 0x17, 0x17, 0x15, 0x00, 0x02, 0x6c, 0x6e, 0x00, 0x01, 0x01, 0x00, 0x08,
            0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x02, 0x65, 0x6f, 0x00, 0x01, 0x01, 0x00, 0x08, 0x01, 0x00, 0x00, 0x00, 0x01,
            0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x02, 0x6f, 0x73, 0x00,
            0x01, 0x01, 0x00, 0x00, 0x01, 0x00, 0x00, 0x08, 0x01, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x80, 0x00, 0x00, 0x00, 0x01, 0x64, 0x00, 0x01, 0x00, 0x80, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x02, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x22, 0x00, 0x81, 0x11, 0x01, 0x7f, 0x00,
        ];

        assert_eq!(decode(&src, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_decode_uncompressed() -> io::Result<()> {
        let src = [
            0x20, // flags = CAT
            0x07, // uncompressed len = 7
            0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73,
        ];

        assert_eq!(decode(&src, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_decode_rle() -> io::Result<()> {
        let src = [
            0x40, // flags = RLE
            0x0d, // uncompressed len = 13
            0x06, 0x06, 0x17, 0x01, 0x07, 0x6f, 0x00, 0x02, 0x01, 0x01, 0x00, 0x00, 0x01, 0x00,
            0x00, 0x0c, 0x02, 0x00, 0x00, 0x08, 0x02, 0x00, 0x00, 0x80, 0x00, 0x00, 0x64, 0x65,
            0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x03, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00,
            0x3a, 0x20, 0x00, 0x00, 0x7c, 0x20, 0x00, 0x00, 0x52, 0x01, 0x00, 0x00, 0x08, 0x04,
            0x00,
        ];

        assert_eq!(decode(&src, 0)?, b"noooooooodles");

        Ok(())
    }

    #[test]
    fn test_decode_bit_packing_with_6_symbols() -> io::Result<()> {
        let src = [
            0x80, // flags = PACK
            0x07, // uncompressed len = 7
            0x06, 0x64, 0x65, 0x6c, 0x6e, 0x6f, 0x73, 0x04, 0x04, 0x05, 0x00, 0x12, 0x43, 0x00,
            0x01, 0x01, 0x01, 0x01, 0x00, 0x0c, 0x02, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x08,
            0x02, 0x00, 0x00, 0x04, 0x02, 0x00,
        ];

        assert_eq!(decode(&src, 0)?, b"noodles");

        Ok(())
    }
}
