use std::io;

use crate::codecs::aac::{
    Model, RangeCoder,
    decode::read_symbol_count,
    rle::{CONTINUE, CONTINUE_CONTEXT, INITIAL_CONTEXT, MODEL_COUNT, MODEL_SYMBOL_COUNT},
};

pub(super) fn decode(src: &mut &[u8], dst: &mut [u8]) -> io::Result<()> {
    let symbol_count = read_symbol_count(src)?;
    let mut models = vec![Model::new(symbol_count); symbol_count.get()];

    let mut rle_models = vec![Model::new(MODEL_SYMBOL_COUNT); MODEL_COUNT];

    let mut coder = RangeCoder::new(src)?;

    let mut iter = dst.iter_mut();
    let mut prev_sym = 0;

    while let Some(d) = iter.next() {
        let sym = models[usize::from(prev_sym)].decode(src, &mut coder)?;

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

        prev_sym = sym;
    }

    Ok(())
}
