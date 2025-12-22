use std::io;

use super::read_symbol_count;
use crate::codecs::aac::{Model, RangeCoder};

pub(super) fn decode(src: &mut &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let symbol_count = read_symbol_count(src)?;
    let max_sym = (symbol_count.get() - 1) as u8;
    let mut model = Model::new(max_sym);

    let mut range_coder = RangeCoder::default();
    range_coder.range_decode_create(src)?;

    for b in dst {
        *b = model.decode(src, &mut range_coder)?;
    }

    Ok(())
}
