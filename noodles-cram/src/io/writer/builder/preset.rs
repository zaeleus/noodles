use crate::container::BlockContentEncoderMap;

#[derive(Default)]
pub enum Preset {
    #[default]
    Medium,
}

impl Preset {
    pub fn records_per_slice(&self) -> usize {
        match self {
            Self::Medium => 10240,
        }
    }

    pub fn block_content_encoder_map(&self) -> BlockContentEncoderMap {
        match self {
            Self::Medium => BlockContentEncoderMap::default(),
        }
    }
}
