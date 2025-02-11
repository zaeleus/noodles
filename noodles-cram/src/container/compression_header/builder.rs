use super::{
    data_series_encodings::DataSeriesEncodings, preservation_map, tag_encodings, CompressionHeader,
};
use crate::io::writer::{Options, Record};

#[derive(Debug, Default)]
pub struct Builder {
    preservation_map_builder: preservation_map::Builder,
    tag_encodings_builder: tag_encodings::Builder,
}

impl Builder {
    pub fn apply_options(&mut self, options: &Options) {
        self.preservation_map_builder.apply_options(options);
    }

    pub fn update(&mut self, record: &Record) {
        self.preservation_map_builder.update(record);
        self.tag_encodings_builder.update(record);
    }

    pub fn build(self) -> CompressionHeader {
        let preservation_map = self.preservation_map_builder.build();
        let data_series_encodings = DataSeriesEncodings::default();
        let tag_encodings = self.tag_encodings_builder.build();
        CompressionHeader::new(preservation_map, data_series_encodings, tag_encodings)
    }
}
