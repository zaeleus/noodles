mod decode;
mod encode;
mod parameter;
mod parameters;

pub use self::{decode::decode, encode::encode};

use super::aac::Model;

struct Models {
    len: Vec<Model>,
    qual: Vec<Model>,
    dup: Model,
    rev: Model,
    sel: Model,
}

impl Models {
    fn new(max_sym: u8, max_sel: u8) -> Models {
        Self {
            len: vec![Model::new(u8::MAX); 4],
            qual: vec![Model::new(max_sym); 1 << 16],
            dup: Model::new(1),
            rev: Model::new(1),
            sel: Model::new(max_sel),
        }
    }
}
