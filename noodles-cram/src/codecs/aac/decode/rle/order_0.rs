use std::io;

use crate::{
    codecs::aac::{Model, RangeCoder},
    io::reader::num::read_u8,
};

pub(super) fn decode(src: &mut &[u8], dst: &mut [u8]) -> io::Result<()> {
    let max_sym = read_u8(src).map(|n| if n == 0 { u8::MAX } else { n - 1 })?;

    let mut model_lit = Model::new(max_sym);
    let mut model_run = vec![Model::new(3); 258];

    let mut range_coder = RangeCoder::default();
    range_coder.range_decode_create(src)?;

    let mut i = 0;

    while i < dst.len() {
        let b = model_lit.decode(src, &mut range_coder)?;
        dst[i] = b;

        let mut part = model_run[usize::from(b)].decode(src, &mut range_coder)?;
        let mut run = usize::from(part);
        let mut rctx = 256;

        while part == 3 {
            part = model_run[rctx].decode(src, &mut range_coder)?;
            rctx = 257;
            run += usize::from(part);
        }

        for j in 1..=run {
            dst[i + j] = b;
        }

        i += run + 1;
    }

    Ok(())
}
