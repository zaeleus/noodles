use std::io;

use super::{CONTINUE, CONTINUE_CONTEXT, INITIAL_CONTEXT, MODEL_COUNT, MODEL_MAX_SYMBOL};
use crate::{
    codecs::aac::{Model, RangeCoder},
    io::reader::num::read_u8,
};

pub(super) fn decode(src: &mut &[u8], dst: &mut [u8]) -> io::Result<()> {
    let max_sym = read_u8(src).map(|n| if n == 0 { u8::MAX } else { n - 1 })?;

    let mut model_lit = vec![Model::new(max_sym); usize::from(max_sym) + 1];
    let mut model_run = vec![Model::new(MODEL_MAX_SYMBOL); MODEL_COUNT];

    let mut range_coder = RangeCoder::default();
    range_coder.range_decode_create(src)?;

    let mut i = 0;
    let mut last = 0;

    while i < dst.len() {
        let b = model_lit[last].decode(src, &mut range_coder)?;
        dst[i] = b;
        last = usize::from(b);

        let mut part = model_run[last].decode(src, &mut range_coder)?;
        let mut run = usize::from(part);
        let mut rctx = INITIAL_CONTEXT;

        while part == CONTINUE {
            part = model_run[rctx].decode(src, &mut range_coder)?;
            rctx = CONTINUE_CONTEXT;
            run += usize::from(part);
        }

        for j in 1..=run {
            dst[i + j] = b;
        }

        i += run + 1;
    }

    Ok(())
}
