pub mod reference;

pub use self::reference::Reference;

#[derive(Debug)]
pub struct Index {
    format: i32,
    reference_sequence_name_index: usize,
    start_position_index: usize,
    end_position_index: usize,
    comment: i32,
    header_line_count: u32,
    names: Vec<String>,
    references: Vec<Reference>,
    unmapped_read_count: Option<u64>,
}

impl Index {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        format: i32,
        reference_sequence_name_index: usize,
        start_position_index: usize,
        end_position_index: usize,
        comment: i32,
        header_line_count: u32,
        names: Vec<String>,
        references: Vec<Reference>,
        unmapped_read_count: Option<u64>,
    ) -> Self {
        Self {
            format,
            reference_sequence_name_index,
            start_position_index,
            end_position_index,
            comment,
            header_line_count,
            names,
            references,
            unmapped_read_count,
        }
    }
}
