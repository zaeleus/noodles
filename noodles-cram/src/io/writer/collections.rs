use std::io::{self, Write};

use super::num::write_int;
use crate::file_definition::Version;

pub(crate) fn write_array<W>(writer: &mut W, version: Version, src: &[u8]) -> io::Result<()>
where
    W: Write,
{
    let len =
        i32::try_from(src.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    write_int(writer, version, len)?;
    writer.write_all(src)?;

    Ok(())
}
