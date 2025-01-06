mod format_version;
pub(crate) mod magic_number;

pub(super) use self::{format_version::read_format_version, magic_number::read_magic_number};
