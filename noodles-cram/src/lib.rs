pub use self::{
    bit_reader::BitReader, container::Container, data_series::DataSeries, reader::Reader,
    record::Record,
};

mod bit_reader;
mod container;
pub mod crai;
mod data_series;
mod huffman;
mod num;
mod rans;
mod reader;
pub mod record;
