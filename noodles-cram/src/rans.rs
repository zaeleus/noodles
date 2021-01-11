use std::{
    convert::TryFrom,
    error, fmt,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::num::read_itf8;

#[derive(Debug, Eq, PartialEq)]
struct TryFromByteError(u8);

impl fmt::Display for TryFromByteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid rANS order: expected 0 or 1, got {}", self.0)
    }
}

impl error::Error for TryFromByteError {}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Order {
    Zero,
    One,
}

impl TryFrom<u8> for Order {
    type Error = TryFromByteError;
    fn try_from(b: u8) -> Result<Self, Self::Error> {
        match b {
            0 => Ok(Self::Zero),
            1 => Ok(Self::One),
            _ => Err(TryFromByteError(b)),
        }
    }
}

pub fn rans_decode<R>(reader: &mut R) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let order = reader.read_u8().and_then(|order| {
        Order::try_from(order).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let _compressed_len = reader.read_u32::<LittleEndian>()?;
    let data_len = reader.read_u32::<LittleEndian>()?;

    let mut buf = vec![0; data_len as usize];

    match order {
        Order::Zero => rans_decode_0(reader, &mut buf)?,
        Order::One => rans_decode_1(reader, &mut buf)?,
    }

    Ok(buf)
}

fn read_frequencies_0<R>(
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

pub fn rans_get_cumulative_freq(r: u32) -> u32 {
    r & 0x0fff
}

pub fn rans_get_symbol_from_freq(cumulative_freqs: &[u32], freq: u32) -> u8 {
    let mut sym = 0;

    while sym < 255 && freq >= cumulative_freqs[(sym + 1) as usize] {
        sym += 1;
    }

    sym
}

pub fn rans_advance_step(r: u32, c: u32, f: u32) -> u32 {
    f * (r >> 12) + (r & 0x0fff) - c
}

pub fn rans_renorm<R>(reader: &mut R, mut r: u32) -> io::Result<u32>
where
    R: Read,
{
    while r < (1 << 23) {
        r = (r << 8) + reader.read_u8().map(u32::from)?;
    }

    Ok(r)
}

pub fn rans_decode_0<R>(reader: &mut R, output: &mut [u8]) -> io::Result<()>
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
        read_frequencies_0(
            reader,
            &mut freqs[sym as usize],
            &mut cumulative_freqs[sym as usize],
        )?;

        if rle > 0 {
            rle -= 1;
            sym += 1;
        } else {
            sym = reader.read_u8()?;

            if sym < 255 && sym == last_sym + 1 {
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

pub fn rans_decode_1<R>(reader: &mut R, output: &mut [u8]) -> io::Result<()>
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rans_decode_with_order_0() -> io::Result<()> {
        let expected = b"noodles";

        let data = vec![
            0x00, 0x25, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x64, 0x82, 0x49, 0x65, 0x00,
            0x82, 0x49, 0x6c, 0x82, 0x49, 0x6e, 0x82, 0x49, 0x6f, 0x00, 0x84, 0x92, 0x73, 0x82,
            0x49, 0x00, 0xe2, 0x06, 0x83, 0x18, 0x74, 0x7b, 0x41, 0x0c, 0x2b, 0xa9, 0x41, 0x0c,
            0x25, 0x31, 0x80, 0x03,
        ];

        let mut reader = &data[..];
        let actual = rans_decode(&mut reader)?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_rans_decode_with_order_1() -> io::Result<()> {
        let expected = b"noodles";

        let data = vec![
            0x01, 0x3b, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x00, 0x64, 0x84, 0x00, 0x6e,
            0x84, 0x00, 0x6f, 0x00, 0x87, 0xff, 0x00, 0x64, 0x6c, 0x8f, 0xff, 0x00, 0x65, 0x00,
            0x73, 0x8f, 0xff, 0x00, 0x6c, 0x65, 0x8f, 0xff, 0x00, 0x6e, 0x6f, 0x8f, 0xff, 0x00,
            0x6f, 0x00, 0x64, 0x87, 0xff, 0x6f, 0x88, 0x00, 0x00, 0x00, 0x00, 0x04, 0x00, 0x02,
            0x02, 0x28, 0x00, 0x01, 0x02, 0x28, 0x00, 0x01, 0x02, 0x60, 0x00, 0x02,
        ];

        let mut reader = &data[..];
        let actual = rans_decode(&mut reader)?;

        assert_eq!(actual, expected);

        Ok(())
    }

    mod order {
        use std::convert::TryFrom;

        use super::super::{Order, TryFromByteError};

        #[test]
        fn test_try_from() {
            assert_eq!(Order::try_from(0), Ok(Order::Zero));
            assert_eq!(Order::try_from(1), Ok(Order::One));
            assert_eq!(Order::try_from(2), Err(TryFromByteError(2)));
        }
    }
}
