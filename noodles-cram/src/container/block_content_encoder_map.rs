//! CRAM container block content-encoder map.

mod builder;

pub use self::builder::Builder;

use std::collections::HashMap;

use crate::{codecs::Encoder, container::block};

/// A CRAM container block content-encoder map.
#[derive(Clone, Debug)]
pub struct BlockContentEncoderMap {
    core_data_encoder: Option<Encoder>,
    data_series_encoders: Vec<Option<Encoder>>,
    tag_values_encoders: HashMap<block::ContentId, Option<Encoder>>,
}

impl BlockContentEncoderMap {
    /// Creates a block content-encoder map builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::container::BlockContentEncoderMap;
    /// let builder = BlockContentEncoderMap::builder();
    /// ````
    pub fn builder() -> Builder {
        Builder::default()
    }

    pub(crate) fn core_data_encoder(&self) -> Option<&Encoder> {
        self.core_data_encoder.as_ref()
    }

    pub(crate) fn data_series_encoders(&self) -> &[Option<Encoder>] {
        &self.data_series_encoders
    }

    pub(crate) fn get_data_series_encoder(
        &self,
        block_content_id: block::ContentId,
    ) -> Option<Option<&Encoder>> {
        let i = (block_content_id as usize) - 1;
        self.data_series_encoders.get(i).map(|e| e.as_ref())
    }

    pub(crate) fn tag_values_encoders(&self) -> &HashMap<block::ContentId, Option<Encoder>> {
        &self.tag_values_encoders
    }

    pub(crate) fn get_tag_values_encoders(
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
        Self::builder().build()
    }
}
