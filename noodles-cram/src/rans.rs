use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::num::read_itf8;

pub fn rans_decode<R>(reader: &mut R) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let order = reader.read_u8()?;

    let _compressed_len = reader.read_u32::<LittleEndian>()?;
    let data_len = reader.read_u32::<LittleEndian>()?;

    let mut buf = vec![0; data_len as usize];

    if order == 0 {
        rans_decode_0(reader, &mut buf)?;
    } else if order == 1 {
        todo!("rans_decode_1");
    } else {
        panic!("invalid order {}", order);
    }

    Ok(buf)
}

fn read_frequencies_0<R>(reader: &mut R, F: &mut [i32], C: &mut [i32]) -> io::Result<()>
where
    R: Read,
{
    let mut sym = reader.read_u8()?;
    let mut last_sym = sym;
    let mut rle = 0;

    loop {
        let f = read_itf8(reader)?;

        F[sym as usize] = f;

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

    C[0] = 0;

    for i in 0..255 {
        C[i + 1] = C[i] + F[i];
    }

    Ok(())
}

pub fn rans_get_cumulative_freq(R: i32) -> i32 {
    R & 0x0fff
}

pub fn rans_get_symbol_from_freq(C: &[i32], f: i32) -> i32 {
    let mut s = 0;

    while f >= C[s + 1] {
        s = s + 1;
    }

    s as i32
}

pub fn rans_advance_step(R: i32, c: i32, f: i32) -> i32 {
    f * (R >> 12) + (R & 0x0fff) - c
}

pub fn rans_renorm<R>(reader: &mut R, mut big_r: i32) -> io::Result<i32>
where
    R: Read,
{
    while big_r < (1 << 23) {
        big_r = (big_r << 8) + reader.read_u8()? as i32;
    }

    Ok(big_r)
}

pub fn rans_decode_0<R>(reader: &mut R, output: &mut [u8]) -> io::Result<()>
where
    R: Read,
{
    let mut F = vec![0i32; 256];
    let mut C = vec![0i32; 256];
    let mut big_r = vec![0i32; 4096];

    read_frequencies_0(reader, &mut F, &mut C)?;

    for j in 0..4 {
        big_r[j] = reader.read_u32::<LittleEndian>()? as i32;
    }

    let mut i = 0;

    while i < output.len() {
        for j in 0..4 {
            if i + j >= output.len() {
                return Ok(());
            }

            let f = rans_get_cumulative_freq(big_r[j]);
            let s = rans_get_symbol_from_freq(&C, f);

            output[i + j] = s as u8;

            big_r[j] = rans_advance_step(big_r[j], C[s as usize], F[s as usize]);
            big_r[j] = rans_renorm(reader, big_r[j])?;
        }

        i += 4;
    }

    Ok(())
}
