use std::{io, num::NonZero};

use crate::{
    codecs::aac::{Model, RangeCoder},
    io::writer::num::write_u8,
};

pub(super) fn encode(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let max_sym = src.iter().max().copied().unwrap_or(0);
    write_u8(dst, max_sym.overflowing_add(1).0)?;

    let symbol_count = NonZero::new(usize::from(max_sym) + 1).unwrap();
    let mut model_lit = vec![Model::new(symbol_count); symbol_count.get()];
    let mut model_run = vec![Model::new(const { NonZero::new(4).unwrap() }); 258];

    let mut range_coder = RangeCoder::default();

    let mut i = 0;
    let mut last = 0;

    while i < src.len() {
        let sym = src[i];
        model_lit[last].encode(dst, &mut range_coder, sym)?;

        let mut run = src[i + 1..].iter().position(|&s| s != sym).unwrap_or(0);
        i += run + 1;

        let mut rctx = usize::from(sym);
        last = usize::from(sym);

        let mut part = run.min(3);
        model_run[rctx].encode(dst, &mut range_coder, part as u8)?;
        rctx = 256;
        run -= part;

        while part == 3 {
            part = run.min(3);
            model_run[rctx].encode(dst, &mut range_coder, part as u8)?;
            rctx = 257;
            run -= part;
        }
    }

    range_coder.range_encode_end(dst)?;

    Ok(())
}
