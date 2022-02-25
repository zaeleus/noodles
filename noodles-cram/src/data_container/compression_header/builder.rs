use crate::{writer::Options, Record};

use super::{
    data_series_encoding_map::DataSeriesEncodingMap, preservation_map, tag_encoding_map,
    CompressionHeader,
};

#[derive(Debug, Default)]
pub struct Builder {
    preservation_map_builder: preservation_map::Builder,
    tag_encoding_map_builder: tag_encoding_map::Builder,
}

impl Builder {
    pub fn apply_options(&mut self, options: &Options) {
        self.preservation_map_builder.apply_options(options);
    }

    pub fn update(&mut self, reference_sequence: &[u8], record: &Record) {
        self.preservation_map_builder
            .update(reference_sequence, record);
        self.tag_encoding_map_builder.update(record);
    }

    pub fn build(self) -> CompressionHeader {
        let preservation_map = self.preservation_map_builder.build();
        let data_series_encoding_map = DataSeriesEncodingMap::default();
        let tag_encoding_map = self.tag_encoding_map_builder.build();
        CompressionHeader::new(preservation_map, data_series_encoding_map, tag_encoding_map)
    }
}
