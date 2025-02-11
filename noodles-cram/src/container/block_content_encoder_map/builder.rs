use std::collections::HashMap;

use super::BlockContentEncoderMap;
use crate::{
    codecs::Encoder,
    container::{
        block,
        compression_header::{data_series_encodings::DataSeries, preservation_map::tag_sets},
    },
};

/// A CRAM container block content-encoder map builder.
#[derive(Debug)]
pub struct Builder {
    core_data_encoder: Option<Encoder>,
    data_series_encoders: Vec<Option<Encoder>>,
    tag_values_encoders: HashMap<block::ContentId, Option<Encoder>>,
}

impl Builder {
    /// Sets the core data encoder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::container::BlockContentEncoderMap;
    /// let builder = BlockContentEncoderMap::builder().set_core_data_encoder(None);
    /// ```
    pub fn set_core_data_encoder(mut self, encoder: Option<Encoder>) -> Self {
        self.core_data_encoder = encoder;
        self
    }

    /// Sets a data series encoder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::container::{
    ///     compression_header::data_series_encodings::DataSeries,
    ///     BlockContentEncoderMap,
    /// };
    ///
    /// let builder = BlockContentEncoderMap::builder()
    ///     .set_data_series_encoder(DataSeries::BamFlags, None);
    /// ```
    pub fn set_data_series_encoder(
        mut self,
        data_series: DataSeries,
        encoder: Option<Encoder>,
    ) -> Self {
        let i = (block::ContentId::from(data_series) as usize) - 1;
        self.data_series_encoders[i] = encoder;
        self
    }

    /// Sets a tag values encoder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::container::{
    ///     compression_header::preservation_map::tag_sets::Key,
    ///     BlockContentEncoderMap,
    /// };
    /// use noodles_sam::alignment::record::data::field::{Tag, Type};
    ///
    /// let key = Key::new(Tag::ALIGNMENT_HIT_COUNT, Type::UInt8);
    /// let builder = BlockContentEncoderMap::builder()
    ///     .set_tag_values_encoder(key, None);
    /// ```
    pub fn set_tag_values_encoder(mut self, key: tag_sets::Key, encoder: Option<Encoder>) -> Self {
        let id = block::ContentId::from(key);
        self.tag_values_encoders.insert(id, encoder);
        self
    }

    /// Builds a block content-encoder map.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::container::BlockContentEncoderMap;
    /// let map = BlockContentEncoderMap::builder().build();
    /// ```
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

        use crate::container::compression_header::data_series_encodings::data_series::STANDARD_DATA_SERIES;

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
