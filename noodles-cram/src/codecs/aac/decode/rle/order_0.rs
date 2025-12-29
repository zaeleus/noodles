use std::io;

use super::{CONTINUE, CONTINUE_CONTEXT, INITIAL_CONTEXT, MODEL_COUNT, MODEL_SYMBOL_COUNT};
use crate::codecs::aac::{Model, RangeCoder, decode::read_symbol_count};

pub(super) fn decode(src: &mut &[u8], dst: &mut [u8]) -> io::Result<()> {
    let symbol_count = read_symbol_count(src)?;
    let mut model = Model::new(symbol_count);

    let mut rle_models = vec![Model::new(MODEL_SYMBOL_COUNT); MODEL_COUNT];

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
