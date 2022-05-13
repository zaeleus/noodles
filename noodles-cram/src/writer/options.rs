#[derive(Clone, Debug)]
pub struct Options {
    pub preserve_read_names: bool,
    pub encode_alignment_start_positions_as_deltas: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            preserve_read_names: true,
            encode_alignment_start_positions_as_deltas: true,
        }
    }
}
