pub(super) mod file_format;
pub(super) mod map;

pub(super) use self::{
    file_format::write_file_format,
    map::{write_map, write_other_map},
};

use std::io::{self, Write};

pub(super) fn write_string<W>(writer: &mut W, s: &str) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(s.as_bytes())
}
