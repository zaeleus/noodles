use std::io::{self, Write};

use crate::file_definition::Version;

pub(super) fn write_format_version<W>(writer: &mut W, version: Version) -> io::Result<()>
where
    W: Write,
{
    let buf = [version.major(), version.minor()];
    writer.write_all(&buf)
}
