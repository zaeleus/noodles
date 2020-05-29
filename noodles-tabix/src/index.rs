mod format;
pub mod reference;

pub use self::{format::Format, reference::Reference};

#[derive(Debug)]
pub struct Index {
    format: Format,
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
        format: Format,
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

    pub fn format(&self) -> Format {
        self.format
    }

    pub fn reference_sequence_name_index(&self) -> usize {
        self.reference_sequence_name_index
    }

    pub fn start_position_index(&self) -> usize {
        self.start_position_index
    }

    pub fn end_position_index(&self) -> usize {
        self.end_position_index
    }

    pub fn comment(&self) -> i32 {
        self.comment
    }

    pub fn header_line_count(&self) -> u32 {
        self.header_line_count
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }

    pub fn references(&self) -> &[Reference] {
        &self.references
    }

    pub fn unmapped_read_count(&self) -> Option<u64> {
        self.unmapped_read_count
    }
}
