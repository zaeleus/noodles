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
        write!(f, "invalid rANS context: expected 0 or 1, got {}", self.0)
    }
}

impl error::Error for TryFromByteError {}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Context {
    Order0,
    Order1,
}

impl TryFrom<u8> for Context {
    type Error = TryFromByteError;

    fn try_from(b: u8) -> Result<Self, Self::Error> {
        match b {
            0 => Ok(Self::Order0),
            1 => Ok(Self::Order1),
            _ => Err(TryFromByteError(b)),
        }
    }
}

pub fn rans_decode<R>(reader: &mut R) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let context = reader.read_u8().and_then(|order| {
        Context::try_from(order).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let _compressed_len = reader.read_u32::<LittleEndian>()?;
    let data_len = reader.read_u32::<LittleEndian>()?;

    let mut buf = vec![0; data_len as usize];

    match context {
        Context::Order0 => {
            rans_decode_0(reader, &mut buf)?;
        }
        Context::Order1 => {
            rans_decode_1(reader, &mut buf)?;
        }
    }

    Ok(buf)
}

fn read_frequencies_0<R>(
    reader: &mut R,
    freqs: &mut [i32],
    cumulative_freqs: &mut [i32],
) -> io::Result<()>
where
    R: Read,
{
    let mut sym = reader.read_u8()?;
    let mut last_sym = sym;
    let mut rle = 0;

    loop {
        let f = read_itf8(reader)?;

        freqs[sym as usize] = f;

        if rle > 0 {
            rle = rle - 1;
            sym = sym + 1;
        } else {
            sym = reader.read_u8()?;

            if sym == last_sym + 1 {
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

pub fn rans_get_cumulative_freq(r: i32) -> i32 {
    r & 0x0fff
}

pub fn rans_get_symbol_from_freq(cumulative_freqs: &[i32], freq: i32) -> i32 {
    let mut sym = 0;

    while freq >= cumulative_freqs[sym + 1] {
        sym += 1;
    }

    sym as i32
}

pub fn rans_advance_step(r: i32, c: i32, f: i32) -> i32 {
    f * (r >> 12) + (r & 0x0fff) - c
}

pub fn rans_renorm<R>(reader: &mut R, mut r: i32) -> io::Result<i32>
where
    R: Read,
{
    while r < (1 << 23) {
        r = (r << 8) + reader.read_u8()? as i32;
    }

    Ok(r)
}

pub fn rans_decode_0<R>(reader: &mut R, output: &mut [u8]) -> io::Result<()>
where
    R: Read,
{
    let mut freqs = vec![0; 256];
    let mut cumulative_freqs = vec![0; 256];
    let mut state = vec![0; 4096];

    read_frequencies_0(reader, &mut freqs, &mut cumulative_freqs)?;

    for j in 0..4 {
        state[j] = reader.read_u32::<LittleEndian>()? as i32;
    }

    let mut i = 0;

    while i < output.len() {
        for j in 0..4 {
            if i + j >= output.len() {
                return Ok(());
            }

            let f = rans_get_cumulative_freq(state[j]);
            let s = rans_get_symbol_from_freq(&cumulative_freqs, f);

            output[i + j] = s as u8;

            state[j] = rans_advance_step(state[j], cumulative_freqs[s as usize], freqs[s as usize]);
            state[j] = rans_renorm(reader, state[j])?;
        }

        i += 4;
    }

    Ok(())
}

fn read_frequencies_1<R>(
    reader: &mut R,
    freqs: &mut Vec<Vec<i32>>,
    cumulative_freqs: &mut Vec<Vec<i32>>,
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
            rle = rle - 1;
            sym = sym + 1;
        } else {
            sym = reader.read_u8()?;

            if sym == last_sym + 1 {
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
    let mut state = vec![0; 4096];
    let mut last_syms = vec![0; 4];

    read_frequencies_1(reader, &mut freqs, &mut cumulative_freqs)?;

    for j in 0..4 {
        state[j] = reader.read_u32::<LittleEndian>()? as i32;
        last_syms[j] = 0;
    }

    let mut i = 0;

    while i < output.len() / 4 {
        for j in 0..4 {
            let f = rans_get_cumulative_freq(state[j]);
            let s = rans_get_symbol_from_freq(&cumulative_freqs[last_syms[j] as usize], f);

            output[i + j * (output.len() / 4)] = s as u8;

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

        output[i + 3 * (output.len() / 4)] = s as u8;

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
    mod context {
        use std::convert::TryFrom;

        use super::super::{Context, TryFromByteError};

        #[test]
        fn test_try_from() {
            assert_eq!(Context::try_from(0), Ok(Context::Order0));
            assert_eq!(Context::try_from(1), Ok(Context::Order1));
            assert_eq!(Context::try_from(2), Err(TryFromByteError(2)));
        }
    }
}
