mod order_0;
mod order_1;

use std::io;

use crate::codecs::aac::Flags;

pub(super) fn decode(src: &mut &[u8], flags: Flags, dst: &mut [u8]) -> io::Result<()> {
    if flags.order() == 0 {
        order_0::decode(src, dst)
    } else {
        order_1::decode(src, dst)
    }
}
