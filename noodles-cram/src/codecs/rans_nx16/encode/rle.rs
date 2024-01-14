use std::io::{self, Write};

use crate::io::writer::num::write_uint7;

pub fn encode(src: &[u8]) -> io::Result<(Vec<u8>, Vec<u8>)> {
    let mut scores = [0; 256];

    for window in src.windows(2) {
        let prev_sym = usize::from(window[0]);
        let curr_sym = usize::from(window[1]);

        if curr_sym == prev_sym {
            scores[curr_sym] += 1;
        } else {
            scores[curr_sym] -= 1;
        }
    }

    let mut n = scores.iter().filter(|&&s| s > 0).count();

    if n == 0 {
        n = 1;
        scores[0] = 1;
    }

    let mut meta = vec![n as u8];

    for (i, &score) in scores.iter().enumerate() {
        if score > 0 {
            let sym = i as u8;
            meta.push(sym);
        }
    }

    let mut buf = vec![0; src.len()];
    let mut end = 0;

    let mut i = 0;

    while i < src.len() {
        buf[end] = src[i];
        end += 1;

        if scores[src[i] as usize] > 0 {
            let mut run = 0;
            let last = src[i];

            while i + run + 1 < src.len() && src[i + run + 1] == last {
                run += 1;
            }

            write_uint7(&mut meta, run as u32)?;

            i += run;
        }

        i += 1;
    }

    buf.truncate(end);

    let mut header = Vec::new();

    let rle_meta_len = meta.len() as u32;
    write_uint7(&mut header, (rle_meta_len << 1) | 1)?;

    let len = buf.len() as u32;
    write_uint7(&mut header, len)?;

    header.write_all(&meta)?;

    Ok((header, buf))
}
