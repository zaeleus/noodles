use std::collections::HashMap;

use super::BlockContentEncoderMap;
use crate::{
    codecs::Encoder,
    container::block,
    data_container::compression_header::{
        data_series_encoding_map::DataSeries, preservation_map::tag_ids_dictionary,
    },
};

#[derive(Debug)]
pub struct Builder {
    core_data_encoder: Option<Encoder>,
    data_series_encoders: Vec<Option<Encoder>>,
    tag_values_encoders: HashMap<block::ContentId, Option<Encoder>>,
}

#[allow(dead_code)]
impl Builder {
    pub fn set_core_data_encoder(mut self, encoder: Option<Encoder>) -> Self {
        self.core_data_encoder = encoder;
        self
    }

    pub fn set_data_series_encoder(
        mut self,
        data_series: DataSeries,
        encoder: Option<Encoder>,
    ) -> Self {
        let i = (i32::from(block::ContentId::from(data_series)) as usize) - 1;
        self.data_series_encoders[i] = encoder;
        self
    }

    pub fn set_tag_values_encoder(
        mut self,
        key: tag_ids_dictionary::Key,
        encoder: Option<Encoder>,
    ) -> Self {
        let id = block::ContentId::from(key);
        self.tag_values_encoders.insert(id, encoder);
        self
    }

    pub fn build(self) -> BlockContentEncoderMap {
        BlockContentEncoderMap {
            core_data_encoder: self.core_data_encoder,
            data_series_encoders: self.data_series_encoders,
            tag_values_encoders: self.tag_values_encoders,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        use flate2::Compression;

        use crate::data_container::compression_header::data_series_encoding_map::data_series::STANDARD_DATA_SERIES;

        let compression_level = Compression::default();

        Self {
            core_data_encoder: Some(Encoder::Gzip(compression_level)),
            data_series_encoders: vec![
                Some(Encoder::Gzip(compression_level));
                STANDARD_DATA_SERIES.len()
            ],
            tag_values_encoders: HashMap::new(),
        }
    }
}
