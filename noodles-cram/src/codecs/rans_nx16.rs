pub mod decode;
mod encode;
mod flags;

pub use self::{decode::rans_decode_nx16, encode::rans_encode_nx16};

use self::flags::Flags;
