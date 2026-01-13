use std::io;

use crate::codecs::aac::{
    Model, RangeCoder,
    encode::{count_symbols, write_symbol_count},
    rle::{CONTINUE, CONTINUE_CONTEXT, INITIAL_CONTEXT, MODEL_COUNT, MODEL_SYMBOL_COUNT},
};

pub(super) fn encode(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let symbol_count = count_symbols(src);
    write_symbol_count(dst, symbol_count)?;

    let mut model_lit = Model::new(symbol_count);
    let mut model_run = vec![Model::new(MODEL_SYMBOL_COUNT); MODEL_COUNT];

    let mut range_coder = RangeCoder::default();

    let mut i = 0;

    while i < src.len() {
        let sym = src[i];
        model_lit.encode(dst, &mut range_coder, sym)?;

        let mut run = src[i + 1..].iter().position(|&s| s != sym).unwrap_or(0);
        i += run + 1;

        let mut rctx = usize::from(sym);

        let mut part = run.min(CONTINUE);
        model_run[rctx].encode(dst, &mut range_coder, part as u8)?;
        rctx = INITIAL_CONTEXT;
        run -= part;

        while part == CONTINUE {
            part = run.min(CONTINUE);
            model_run[rctx].encode(dst, &mut range_coder, part as u8)?;
            rctx = CONTINUE_CONTEXT;
            run -= part;
        }
    }

    range_coder.range_encode_end(dst)?;

    Ok(())
}
