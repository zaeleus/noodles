pub use self::{
    bit_reader::BitReader, container::Container, data_series::DataSeries, encoding::Encoding,
    reader::Reader, record::Record,
};

mod bit_reader;
mod container;
pub mod crai;
mod data_series;
mod encoding;
mod huffman;
mod num;
mod rans;
mod reader;
pub mod record;

use std::collections::HashMap;

use crate::num::Itf8;

pub type TagEncodingMap = HashMap<Itf8, Encoding>;
