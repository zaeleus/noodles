pub use self::reader::block::read_block;
pub use self::{
    bit_reader::BitReader, block::Block, container::Container, data_series::DataSeries,
    encoding::Encoding, reader::Reader, record::Record,
};

mod bit_reader;
mod block;
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

pub type DataSeriesEncodingMap = HashMap<DataSeries, Encoding>;
pub type TagEncodingMap = HashMap<Itf8, Encoding>;
