use std::io;

use super::{CONTINUE, CONTINUE_CONTEXT, INITIAL_CONTEXT, MODEL_COUNT, MODEL_MAX_SYMBOL};
use crate::{
    codecs::aac::{Model, RangeCoder},
    io::reader::num::read_u8,
};

pub(super) fn decode(src: &mut &[u8], dst: &mut [u8]) -> io::Result<()> {
    let max_sym = read_u8(src).map(|n| if n == 0 { u8::MAX } else { n - 1 })?;

    let mut model = Model::new(max_sym);
    let mut rle_models = vec![Model::new(MODEL_MAX_SYMBOL); MODEL_COUNT];

    let mut coder = RangeCoder::default();
    coder.range_decode_create(src)?;

    let mut iter = dst.iter_mut();

    while let Some(d) = iter.next() {
        let sym = model.decode(src, &mut coder)?;

        *d = sym;

        let mut rle_model = &mut rle_models[usize::from(sym)];

        let mut n = rle_model.decode(src, &mut coder)?;
        let mut len = usize::from(n);

        rle_model = &mut rle_models[INITIAL_CONTEXT];

        while n == CONTINUE {
            n = rle_model.decode(src, &mut coder)?;
            len += usize::from(n);
            rle_model = &mut rle_models[CONTINUE_CONTEXT];
        }

        for e in iter.by_ref().take(len) {
            *e = sym;
        }
    }

    Ok(())
}
