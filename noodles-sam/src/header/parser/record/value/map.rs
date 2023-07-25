mod field;
pub(crate) mod header;
pub(crate) mod program;
pub(crate) mod read_group;
pub(crate) mod reference_sequence;

pub(crate) use self::{
    header::parse_header, program::parse_program, read_group::parse_read_group,
    reference_sequence::parse_reference_sequence,
};
