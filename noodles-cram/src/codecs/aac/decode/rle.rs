mod order_0;
mod order_1;

use std::{io, num::NonZero};

use crate::codecs::aac::Flags;

const INITIAL_CONTEXT: usize = 256;
const CONTINUE_CONTEXT: usize = INITIAL_CONTEXT + 1;

const CONTINUE: u8 = 3;

const MODEL_SYMBOL_COUNT: NonZero<usize> = NonZero::new(4).unwrap();
const MODEL_COUNT: usize = CONTINUE_CONTEXT + 1;

pub(super) fn decode(src: &mut &[u8], flags: Flags, dst: &mut [u8]) -> io::Result<()> {
    if flags.order() == 0 {
        order_0::decode(src, dst)
    } else {
        order_1::decode(src, dst)
    }
}
