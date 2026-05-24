use std::array;

use crate::{codecs::Encoder, container::BlockContentEncoderMap};

#[derive(Default)]
pub enum Preset {
    #[default]
    Medium,
    #[expect(dead_code)]
    Slow,
}

impl Preset {
    pub fn records_per_slice(&self) -> usize {
        match self {
            Self::Medium => 10240,
            Self::Slow => 25600,
        }
    }

    pub fn block_content_encoder_map(&self) -> BlockContentEncoderMap {
        match self {
            Self::Medium => BlockContentEncoderMap::default(),
            Self::Slow => {
                let encoder = Some(Encoder::Bzip2(Default::default()));

                BlockContentEncoderMap {
                    core_data_encoder: encoder.clone(),
                    data_series_encoders: array::from_fn(|_| encoder.clone()),
                    tag_values_encoders: Default::default(),
                    default_encoder: encoder,
                }
            }
        }
    }
}
