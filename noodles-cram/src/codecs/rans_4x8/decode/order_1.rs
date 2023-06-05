use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use super::{order_0, rans_advance_step, rans_get_cumulative_freq, rans_renorm};

pub fn decode<R>(reader: &mut R, dst: &mut [u8]) -> io::Result<()>
where
    R: Read,
{
    let mut freqs = [[0; 256]; 256];
    let mut cumulative_freqs = [[0; 256]; 256];

    read_frequencies_1(reader, &mut freqs, &mut cumulative_freqs)?;

    let cumulative_freqs_symbols_tables = build_cumulative_freqs_symbols_table_1(&cumulative_freqs);

    let mut states = [0; 4];
    reader.read_u32_into::<LittleEndian>(&mut states)?;

    let state_count = states.len();
    let chunk_size = dst.len() / state_count;
    let (left, right) = dst.split_at_mut(2 * chunk_size);
    let (chunk_0, chunk_1) = left.split_at_mut(chunk_size);
    let (chunk_2, chunk_3) = right.split_at_mut(chunk_size);
    let mut chunks = [
        (states[0], 0, chunk_0),
        (states[1], 0, chunk_1),
        (states[2], 0, chunk_2),
        (states[3], 0, chunk_3),
    ];

    for i in 0..chunk_size {
        for (r, last_sym, chunk) in &mut chunks {
            let f = rans_get_cumulative_freq(*r);
            let s = cumulative_freqs_symbols_tables[usize::from(*last_sym)][f as usize];

            chunk[i] = s;

            *r = rans_advance_step(
                *r,
                cumulative_freqs[usize::from(*last_sym)][usize::from(s)],
                freqs[usize::from(*last_sym)][usize::from(s)],
            );
            *r = rans_renorm(reader, *r)?;

            *last_sym = s;
        }
    }

    let (mut r, mut last_sym, chunk) = &mut chunks[3];
    let remainder = &mut chunk[chunk_size..];

    for d in remainder {
        let f = rans_get_cumulative_freq(r);
        let s = cumulative_freqs_symbols_tables[usize::from(last_sym)][f as usize];

        *d = s;

        r = rans_advance_step(
            r,
            cumulative_freqs[usize::from(last_sym)][usize::from(s)],
            freqs[usize::from(last_sym)][usize::from(s)],
        );
        r = rans_renorm(reader, r)?;

        last_sym = s;
    }

    Ok(())
}

fn read_frequencies_1<R>(
    reader: &mut R,
    freqs: &mut [[u32; 256]; 256],
    cumulative_freqs: &mut [[u32; 256]; 256],
) -> io::Result<()>
where
    R: Read,
{
    let mut sym = reader.read_u8()?;
    let mut last_sym = sym;
    let mut rle = 0;

    loop {
        order_0::read_frequencies_0(
            reader,
            &mut freqs[usize::from(sym)],
            &mut cumulative_freqs[usize::from(sym)],
        )?;

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

    Ok(())
}

pub fn build_cumulative_freqs_symbols_table_1(
    cumulative_freqs: &[[u32; 256]; 256],
) -> Box<[[u8; 4096]; 256]> {
    let mut tables = Box::new([[0; 4096]; 256]);

    for (table, cumulative_freqs) in tables.iter_mut().zip(cumulative_freqs) {
        *table = order_0::build_cumulative_freqs_symbols_table_0(cumulative_freqs);
    }

    tables
}
