use crate::{data_container::BlockContentEncoderMap, file_definition::Version};

#[derive(Clone, Debug)]
pub struct Options {
    pub preserve_read_names: bool,
    pub encode_alignment_start_positions_as_deltas: bool,
    pub version: Version,
    pub block_content_encoder_map: BlockContentEncoderMap,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            preserve_read_names: true,
            encode_alignment_start_positions_as_deltas: true,
            version: Version::default(),
            block_content_encoder_map: BlockContentEncoderMap::default(),
        }
    }
}
