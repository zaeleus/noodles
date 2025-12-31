use std::num::NonZero;

use crate::codecs::aac::Model;

// § 6.1 "FQZComp Models" (2023-03-15).
const QUALITY_MODEL_COUNT: usize = 1 << 16;

const LENGTH_MODEL_SYMBOL_COUNT: NonZero<usize> = NonZero::new(256).unwrap();
const LENGTH_MODEL_COUNT: usize = 4;

const BOOL_MODEL_SYMBOL_COUNT: NonZero<usize> = NonZero::new(2).unwrap();
const REVERSE_MODEL_SYMBOL_COUNT: NonZero<usize> = BOOL_MODEL_SYMBOL_COUNT;
const DUPLICATE_MODEL_SYMBOL_COUNT: NonZero<usize> = BOOL_MODEL_SYMBOL_COUNT;

pub struct Models {
    pub qual: Vec<Model>,
    pub len: Vec<Model>,
    pub rev: Model,
    pub dup: Model,
    pub sel: Option<Model>,
}

impl Models {
    pub fn new(max_symbol_count: NonZero<usize>, selector_count: Option<NonZero<usize>>) -> Models {
        Self {
            qual: vec![Model::new(max_symbol_count); QUALITY_MODEL_COUNT],
            len: vec![Model::new(LENGTH_MODEL_SYMBOL_COUNT); LENGTH_MODEL_COUNT],
            rev: Model::new(REVERSE_MODEL_SYMBOL_COUNT),
            dup: Model::new(DUPLICATE_MODEL_SYMBOL_COUNT),
            sel: selector_count.map(Model::new),
        }
    }
}
