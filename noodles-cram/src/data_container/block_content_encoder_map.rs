use std::collections::HashMap;

use crate::{codecs::Encoder, container::block};

#[derive(Clone, Debug)]
pub struct BlockContentEncoderMap {
    core_data_encoder: Option<Encoder>,
    data_series_encoders: Vec<Option<Encoder>>,
    tag_values_encoders: HashMap<block::ContentId, Option<Encoder>>,
}

impl BlockContentEncoderMap {
    pub fn core_data_encoder(&self) -> Option<&Encoder> {
        self.core_data_encoder.as_ref()
    }

    pub fn get_data_series_encoder(
        &self,
        block_content_id: block::ContentId,
    ) -> Option<Option<&Encoder>> {
        let i = (i32::from(block_content_id) as usize) - 1;
        self.data_series_encoders.get(i).map(|e| e.as_ref())
    }

    pub fn get_tag_values_encoders(
        &self,
        block_content_id: block::ContentId,
    ) -> Option<Option<&Encoder>> {
        self.tag_values_encoders
            .get(&block_content_id)
            .map(|e| e.as_ref())
    }
}

impl Default for BlockContentEncoderMap {
    fn default() -> Self {
        use super::compression_header::data_series_encoding_map::data_series::STANDARD_DATA_SERIES;

        Self {
            core_data_encoder: Some(Encoder::Gzip),
            data_series_encoders: vec![Some(Encoder::Gzip); STANDARD_DATA_SERIES.len()],
            tag_values_encoders: HashMap::new(),
        }
    }
}
