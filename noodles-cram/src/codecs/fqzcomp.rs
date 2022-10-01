mod decode;
mod encode;
mod parameter;
mod parameters;

pub use self::{decode::decode, encode::encode};

use self::parameters::Parameters;
use super::aac::{Model, RangeCoder};

struct Models {
    len: Vec<Model>,
    qual: Vec<Model>,
    dup: Model,
    rev: Model,
    sel: Model,
}

fn fqz_create_models(parameters: &Parameters) -> (RangeCoder, Models) {
    let models = Models {
        len: vec![Model::new(u8::MAX); 4],
        qual: vec![Model::new(parameters.max_sym); 1 << 16],
        dup: Model::new(1),
        rev: Model::new(1),
        sel: Model::new(parameters.max_sel),
    };

    (RangeCoder::default(), models)
}
