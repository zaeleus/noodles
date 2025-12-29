use std::num::NonZero;

use crate::codecs::aac::Model;

pub struct Models {
    pub len: Vec<Model>,
    pub qual: Vec<Model>,
    pub dup: Model,
    pub rev: Model,
    pub sel: Option<Model>,
}

impl Models {
    pub fn new(max_symbol_count: NonZero<usize>, selector_count: Option<NonZero<usize>>) -> Models {
        Self {
            len: vec![Model::new(const { NonZero::new(256).unwrap() }); 4],
            qual: vec![Model::new(max_symbol_count); 1 << 16],
            dup: Model::new(const { NonZero::new(2).unwrap() }),
            rev: Model::new(const { NonZero::new(2).unwrap() }),
            sel: selector_count.map(Model::new),
        }
    }
}
