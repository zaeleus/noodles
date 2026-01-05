use std::io;

use super::read_symbol_count;
use crate::codecs::aac::{Model, RangeCoder};

pub(super) fn decode(src: &mut &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let symbol_count = read_symbol_count(src)?;
    let mut model = Model::new(symbol_count);

    let mut range_coder = RangeCoder::new(src)?;

    for b in dst {
        *b = model.decode(src, &mut range_coder)?;
    }

    Ok(())
}
