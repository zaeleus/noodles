mod builder;
pub mod data_series_encoding_map;
pub mod encoding;
pub mod preservation_map;
mod tag_encoding_map;

pub use self::{
    builder::Builder,
    data_series_encoding_map::DataSeriesEncodingMap,
    encoding::Encoding,
    preservation_map::{PreservationMap, SubstitutionMatrix, TagIdsDictionary},
    tag_encoding_map::TagEncodingMap,
};

use std::{convert::TryFrom, io};

use crate::reader::compression_header::read_compression_header;

use super::Block;

/// A CRAM data container compression header.
///
/// The compression header has three maps with information about how the data is compressed: a
/// preservation map, a data series encoding map, and a tag encoding map.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompressionHeader {
    preservation_map: PreservationMap,
    data_series_encoding_map: DataSeriesEncodingMap,
    tag_encoding_map: TagEncodingMap,
}

impl CompressionHeader {
    pub fn builder() -> Builder {
        Builder::default()
    }

    pub fn new(
        preservation_map: PreservationMap,
        data_series_encoding_map: DataSeriesEncodingMap,
        tag_encoding_map: TagEncodingMap,
    ) -> Self {
        Self {
            preservation_map,
            data_series_encoding_map,
            tag_encoding_map,
        }
    }

    pub fn preservation_map(&self) -> &PreservationMap {
        &self.preservation_map
    }

    pub fn data_series_encoding_map(&self) -> &DataSeriesEncodingMap {
        &self.data_series_encoding_map
    }

    pub fn tag_encoding_map(&self) -> &TagEncodingMap {
        &self.tag_encoding_map
    }
}

impl TryFrom<&Block> for CompressionHeader {
    type Error = io::Error;

    fn try_from(block: &Block) -> Result<Self, Self::Error> {
        let data = block.decompressed_data()?;
        let mut reader = &data[..];
        read_compression_header(&mut reader)
    }
}
