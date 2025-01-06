mod file_id;
mod format_version;
mod magic_number;

pub(super) use self::{
    file_id::read_file_id, format_version::read_format_version, magic_number::read_magic_number,
};
