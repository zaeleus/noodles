#![allow(dead_code)]

mod flags;

use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use self::flags::Flags;
use crate::reader::num::read_uint7;

pub fn rans_decode_nx16<R>(reader: &mut R, mut len: usize) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let flags = reader.read_u8().map(Flags::from)?;

    if !flags.contains(Flags::NO_SIZE) {
        len = read_uint7(reader).and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;
    }

    if flags.contains(Flags::STRIPE) {
        todo!("rans_decode_stripe");
    }

    let n = if flags.contains(Flags::N32) { 32 } else { 4 };

    if flags.contains(Flags::PACK) {
        todo!("decode_pack_meta");
    }

    if flags.contains(Flags::RLE) {
        todo!("decode_rle_meta");
    }

    let mut data = vec![0; len];

    if flags.contains(Flags::CAT) {
        reader.read_exact(&mut data)?;
    } else if flags.contains(Flags::ORDER) {
        todo!("rans_decode_nx16_1");
    } else {
        rans_decode_nx16_0(reader, &mut data, n)?;
    };

    if flags.contains(Flags::RLE) {
        todo!("decode_rle");
    }

    if flags.contains(Flags::PACK) {
        todo!("decode_pack");
    }

    Ok(data)
}

fn read_alphabet<R>(reader: &mut R) -> io::Result<[bool; 256]>
where
    R: Read,
{
    let mut alphabet = [false; 256];

    let mut sym = reader.read_u8()?;
    let mut last_sym = sym;
    let mut rle = 0;

    loop {
        alphabet[usize::from(sym)] = true;

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

    Ok(alphabet)
}

fn normalize_frequencies_nx16_0(freqs: &mut [u32], bits: u32) {
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

fn read_frequencies_nx16_0<R>(
    reader: &mut R,
    freqs: &mut [u32],
    cumulative_freqs: &mut [u32],
) -> io::Result<()>
where
    R: Read,
{
    let alphabet = read_alphabet(reader)?;

    for i in 0..alphabet.len() {
        if alphabet[i] {
            freqs[i] = read_uint7(reader)?;
        }
    }

    normalize_frequencies_nx16_0(freqs, 12);

    cumulative_freqs[0] = 0;

    for i in 0..255 {
        cumulative_freqs[i + 1] = cumulative_freqs[i] + freqs[i];
    }

    Ok(())
}

fn rans_get_cumulative_freq_nx16(r: u32, bits: u32) -> u32 {
    r & ((1 << bits) - 1)
}

fn rans_get_symbol_from_freq(cumulative_freqs: &[u32], freq: u32) -> u8 {
    let mut sym = 0;

    while sym < 255 && freq >= cumulative_freqs[(sym + 1) as usize] {
        sym += 1;
    }

    sym
}

fn rans_advance_step_nx16(r: u32, c: u32, f: u32, bits: u32) -> u32 {
    f * (r >> bits) + (r & ((1 << bits) - 1)) - c
}

pub fn rans_renorm_nx16<R>(reader: &mut R, mut r: u32) -> io::Result<u32>
where
    R: Read,
{
    if r < (1 << 15) {
        r = (r << 16) + reader.read_u16::<LittleEndian>().map(u32::from)?;
    }

    Ok(r)
}

fn rans_decode_nx16_0<R>(reader: &mut R, output: &mut [u8], n: u32) -> io::Result<()>
where
    R: Read,
{
    let mut freqs = [0; 256];
    let mut cumulative_freqs = [0; 256];

    read_frequencies_nx16_0(reader, &mut freqs, &mut cumulative_freqs)?;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rans_decode_nx16_order_0() -> io::Result<()> {
        let data = [
            0x00, // flags = {empty}
            0x07, // uncompressed len = 7
            0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x01, 0x01, 0x01, 0x01, 0x03,
            0x01, 0x00, 0x26, 0x20, 0x00, 0x00, 0xb8, 0x0a, 0x00, 0x00, 0xd8, 0x0a, 0x00, 0x00,
            0x00, 0x04, 0x00,
        ];
        let mut reader = &data[..];

        assert_eq!(rans_decode_nx16(&mut reader, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_rans_decode_nx16_uncompressed() -> io::Result<()> {
        let data = [
            0x20, // flags = CAT
            0x07, // uncompressed len = 7
            0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73,
        ];

        let mut reader = &data[..];
        assert_eq!(rans_decode_nx16(&mut reader, 0)?, b"noodles");

        Ok(())
    }
}
