//! CRAM data container compress header.

mod builder;
pub mod data_series_encoding_map;
pub(crate) mod encoding;
pub mod preservation_map;
mod tag_encoding_map;

pub(crate) use self::{
    builder::Builder,
    data_series_encoding_map::DataSeriesEncodingMap,
    encoding::Encoding,
    preservation_map::{PreservationMap, SubstitutionMatrix, TagSets},
    tag_encoding_map::TagEncodingMap,
};

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
    pub(crate) fn builder() -> Builder {
        Builder::default()
    }

    pub(crate) fn new(
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

    pub(crate) fn preservation_map(&self) -> &PreservationMap {
        &self.preservation_map
    }

    pub(crate) fn data_series_encoding_map(&self) -> &DataSeriesEncodingMap {
        &self.data_series_encoding_map
    }

    pub(crate) fn tag_encoding_map(&self) -> &TagEncodingMap {
        &self.tag_encoding_map
    }
}
