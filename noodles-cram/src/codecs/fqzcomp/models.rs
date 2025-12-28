use crate::codecs::aac::Model;

pub struct Models {
    pub len: Vec<Model>,
    pub qual: Vec<Model>,
    pub dup: Model,
    pub rev: Model,
    pub sel: Option<Model>,
}

impl Models {
    pub fn new(max_sym: u8, max_sel: Option<u8>) -> Models {
        Self {
            len: vec![Model::new(u8::MAX); 4],
            qual: vec![Model::new(max_sym); 1 << 16],
            dup: Model::new(1),
            rev: Model::new(1),
            sel: max_sel.map(Model::new),
        }
    }
}
