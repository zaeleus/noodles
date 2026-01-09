use std::num::NonZero;

pub const INITIAL_CONTEXT: usize = 256;
pub const CONTINUE_CONTEXT: usize = INITIAL_CONTEXT + 1;

pub const CONTINUE: u8 = 3;

pub const MODEL_SYMBOL_COUNT: NonZero<usize> = NonZero::new(4).unwrap();
pub const MODEL_COUNT: usize = CONTINUE_CONTEXT + 1;
