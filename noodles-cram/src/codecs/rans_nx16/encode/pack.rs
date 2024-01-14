use std::io;

use byteorder::WriteBytesExt;

use crate::io::writer::num::write_uint7;

pub fn encode(src: &[u8]) -> io::Result<(Vec<u8>, Vec<u8>)> {
    let mut frequencies = [0; 256];

    for &b in src {
        let sym = usize::from(b);
        frequencies[sym] += 1;
    }

    let mut lut = [0; 256];
    let mut n = 0;

    for (sym, &f) in frequencies.iter().enumerate() {
        if f > 0 {
            lut[sym] = n;
            n += 1;
        }
    }

    let buf = if n <= 1 {
        Vec::new()
    } else if n <= 2 {
        let len = (src.len() / 8) + 1;
        let mut dst = vec![0; len];

        for (d, chunk) in dst.iter_mut().zip(src.chunks(8)) {
            for (shift, &s) in chunk.iter().enumerate() {
                let sym = usize::from(s);
                let value = lut[sym];
                *d |= value << shift;
            }
        }

        dst
    } else if n <= 4 {
        let len = (src.len() / 4) + 1;
        let mut dst = vec![0; len];

        for (d, chunk) in dst.iter_mut().zip(src.chunks(4)) {
            for (shift, &s) in chunk.iter().enumerate() {
                let sym = usize::from(s);
                let value = lut[sym];
                *d |= value << (shift * 2);
            }
        }

        dst
    } else if n <= 16 {
        let len = (src.len() / 2) + 1;
        let mut dst = vec![0; len];

        for (d, chunk) in dst.iter_mut().zip(src.chunks(2)) {
            for (shift, &s) in chunk.iter().enumerate() {
                let sym = usize::from(s);
                let value = lut[sym];
                *d |= value << (shift * 4);
            }
        }

        dst
    } else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "unique symbols > 16",
        ));
    };

    let mut header = Vec::new();
    header.write_u8(n)?;

    for (sym, &f) in frequencies.iter().enumerate() {
        if f > 0 {
            let b = sym as u8;
            header.write_u8(b)?;
        }
    }

    let len = buf.len() as u32;
    write_uint7(&mut header, len)?;

    Ok((header, buf))
}
