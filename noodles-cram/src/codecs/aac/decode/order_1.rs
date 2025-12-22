use std::io;

use super::read_symbol_count;
use crate::codecs::aac::{Model, RangeCoder};

pub(super) fn decode(src: &mut &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let symbol_count = read_symbol_count(src)?;
    let max_sym = (symbol_count.get() - 1) as u8;
    let mut models = vec![Model::new(max_sym); usize::from(max_sym) + 1];

    let mut range_coder = RangeCoder::default();
    range_coder.range_decode_create(src)?;

    let mut last = 0;

    for b in dst {
        *b = models[last].decode(src, &mut range_coder)?;
        last = usize::from(*b);
    }

    Ok(())
}
