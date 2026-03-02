use crate::{container::BlockContentEncoderMap, file_definition::Version};

pub(crate) const DEFAULT_RECORDS_PER_SLICE: usize = 10240;
pub(crate) const DEFAULT_SLICES_PER_CONTAINER: usize = 1;

#[derive(Clone, Debug)]
pub struct Options {
    pub preserve_read_names: bool,
    pub encode_alignment_start_positions_as_deltas: bool,
    pub version: Version,
    pub block_content_encoder_map: BlockContentEncoderMap,
    pub records_per_slice: usize,
    pub slices_per_container: usize,
    pub embed_reference_sequences: bool,
    pub strip_md_nm: bool,
    pub reference_required: bool,
    /// CRAM 4.0 quality score orientation: `true` = alignment orientation (QO=1),
    /// `false` = original/sequencing orientation (QO=0, requires reversal for
    /// reverse-strand reads). Ignored for CRAM 2.x/3.x.
    pub qs_seq_orient: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            preserve_read_names: true,
            encode_alignment_start_positions_as_deltas: true,
            version: Version::default(),
            block_content_encoder_map: BlockContentEncoderMap::default(),
            records_per_slice: DEFAULT_RECORDS_PER_SLICE,
            slices_per_container: DEFAULT_SLICES_PER_CONTAINER,
            embed_reference_sequences: false,
            strip_md_nm: false,
            reference_required: true,
            qs_seq_orient: true,
        }
    }
}
