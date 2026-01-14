use std::io;

use crate::codecs::aac::{
    Model, RangeCoder,
    encode::{count_symbols, write_symbol_count},
    rle::{CONTINUE, CONTINUE_CONTEXT, INITIAL_CONTEXT, MODEL_COUNT, MODEL_SYMBOL_COUNT},
};

pub(super) fn encode(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let symbol_count = count_symbols(src);
    write_symbol_count(dst, symbol_count)?;

    let mut model = Model::new(symbol_count);
    let mut rle_models = vec![Model::new(MODEL_SYMBOL_COUNT); MODEL_COUNT];

    let mut coder = RangeCoder::default();

    let mut iter = src.iter().peekable();

    while let Some(&sym) = iter.next() {
        model.encode(dst, &mut coder, sym)?;

        let mut len = 0;

        while let Some(&&s) = iter.peek()
            && s == sym
        {
            len += 1;
            iter.next();
        }

        let mut rle_model = &mut rle_models[usize::from(sym)];

        let mut n = len.min(CONTINUE);
        // SAFETY: `n <= CONTINUE`.
        rle_model.encode(dst, &mut coder, n as u8)?;
        len -= n;
        rle_model = &mut rle_models[INITIAL_CONTEXT];

        while n == CONTINUE {
            n = len.min(CONTINUE);
            // SAFETY: `n <= CONTINUE`.
            rle_model.encode(dst, &mut coder, n as u8)?;
            len -= n;
            rle_model = &mut rle_models[CONTINUE_CONTEXT];
        }
    }

    coder.range_encode_end(dst)?;

    Ok(())
}
