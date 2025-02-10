use std::io::{self, Write};

use super::num::write_itf8;

pub(crate) fn write_array<W>(writer: &mut W, src: &[u8]) -> io::Result<()>
where
    W: Write,
{
    let len =
        i32::try_from(src.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    write_itf8(writer, len)?;
    writer.write_all(src)?;

    Ok(())
}
