//! Async BGZF.

mod block_codec;
pub mod fs;
pub mod io;

use self::block_codec::BlockCodec;
