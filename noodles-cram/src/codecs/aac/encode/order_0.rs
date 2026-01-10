use std::{io, num::NonZero};

use crate::{
    codecs::aac::{Model, RangeCoder},
    io::writer::num::write_u8,
};

pub(super) fn encode(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let max_sym = src.iter().max().copied().unwrap_or(0);
    write_u8(dst, max_sym.overflowing_add(1).0)?;

    let symbol_count = NonZero::new(usize::from(max_sym) + 1).unwrap();
    let mut model = Model::new(symbol_count);
    let mut range_coder = RangeCoder::default();

    for &sym in src {
        model.encode(dst, &mut range_coder, sym)?;
    }

    range_coder.range_encode_end(dst)?;

    Ok(())
}
