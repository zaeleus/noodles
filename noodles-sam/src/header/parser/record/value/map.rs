mod field;
pub(crate) mod header;
pub(crate) mod reference_sequence;

pub(crate) use self::{header::parse_header, reference_sequence::parse_reference_sequence};
