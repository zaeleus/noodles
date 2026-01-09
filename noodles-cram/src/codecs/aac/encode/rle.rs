mod order_0;
mod order_1;

use std::io;

use crate::codecs::aac::Flags;

pub(super) fn encode(src: &[u8], flags: Flags, dst: &mut Vec<u8>) -> io::Result<()> {
    if flags.order() == 0 {
        order_0::encode(src, dst)
    } else {
        order_1::encode(src, dst)
    }
}
