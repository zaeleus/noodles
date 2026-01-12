use std::io;

use super::{count_symbols, write_symbol_count};
use crate::codecs::aac::{Model, RangeCoder};

pub(super) fn encode(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let symbol_count = count_symbols(src);
    write_symbol_count(dst, symbol_count)?;

    let mut model = Model::new(symbol_count);
    let mut range_coder = RangeCoder::default();

    for &sym in src {
        model.encode(dst, &mut range_coder, sym)?;
    }

    range_coder.range_encode_end(dst)?;

    Ok(())
}
