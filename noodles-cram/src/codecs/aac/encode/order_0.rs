use std::{io, num::NonZero};

use super::write_symbol_count;
use crate::codecs::aac::{Model, RangeCoder};

pub(super) fn encode(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let max_sym = src.iter().max().copied().unwrap_or(0);
    let symbol_count = NonZero::new(usize::from(max_sym) + 1).unwrap();
    write_symbol_count(dst, symbol_count)?;

    let mut model = Model::new(symbol_count);
    let mut range_coder = RangeCoder::default();

    for &sym in src {
        model.encode(dst, &mut range_coder, sym)?;
    }

    range_coder.range_encode_end(dst)?;

    Ok(())
}
