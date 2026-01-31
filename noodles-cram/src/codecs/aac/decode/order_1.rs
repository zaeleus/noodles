use std::io;

use super::read_symbol_count;
use crate::codecs::aac::{Model, RangeCoder};

pub(super) fn decode(src: &mut &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let symbol_count = read_symbol_count(src)?;
    let mut models = vec![Model::new(symbol_count); symbol_count.get()];

    let mut range_coder = RangeCoder::new(src)?;

    let mut prev_sym = 0;

    for d in dst {
        *d = models[usize::from(prev_sym)].decode(src, &mut range_coder)?;
        prev_sym = *d;
    }

    Ok(())
}
