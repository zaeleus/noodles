pub use self::reader::block::read_block;
pub use self::{
    bit_reader::BitReader, block::Block, compression_header::CompressionHeader,
    container::Container, data_series::DataSeries, encoding::Encoding, feature::Feature,
    preservation_map::PreservationMap, reader::Reader, record::Record, slice::Slice, tag::Tag,
};

mod bit_reader;
mod block;
mod compression_header;
mod container;
mod data_series;
mod encoding;
mod feature;
mod huffman;
mod num;
mod preservation_map;
mod rans;
mod reader;
mod record;
mod slice;
mod tag;

use std::collections::HashMap;

use crate::num::Itf8;

pub type DataSeriesEncodingMap = HashMap<DataSeries, Encoding>;
pub type TagEncodingMap = HashMap<Itf8, Encoding>;
