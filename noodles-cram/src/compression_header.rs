use crate::{DataSeriesEncodingMap, PreservationMap, TagEncodings};

#[derive(Debug, Default)]
pub struct CompressionHeader {
    preservation_map: PreservationMap,
    data_series_encoding_map: DataSeriesEncodingMap,
    tag_encodings: TagEncodings,
}

impl CompressionHeader {
    pub fn new(
        preservation_map: PreservationMap,
        data_series_encoding_map: DataSeriesEncodingMap,
        tag_encodings: TagEncodings,
    ) -> Self {
        Self {
            preservation_map,
            data_series_encoding_map,
            tag_encodings,
        }
    }
}
