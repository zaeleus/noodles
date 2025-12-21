use std::io;

use crate::{
    codecs::aac::{Model, RangeCoder},
    io::reader::num::read_u8,
};

pub(super) fn decode(src: &mut &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let max_sym = read_u8(src).map(|n| n.overflowing_sub(1).0)?;

    let mut model = Model::new(max_sym);

    let mut range_coder = RangeCoder::default();
    range_coder.range_decode_create(src)?;

    for b in dst {
        *b = model.decode(src, &mut range_coder)?;
    }

    Ok(())
}
