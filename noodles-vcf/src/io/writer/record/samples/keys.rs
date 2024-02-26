use std::io::{self, Write};

use crate::variant::record_buf::samples::Keys;

pub(super) fn write_keys<W>(writer: &mut W, keys: &Keys) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b":";

    for (i, key) in keys.iter().enumerate() {
        if i > 0 {
            writer.write_all(DELIMITER)?;
        }

        writer.write_all(key.as_bytes())?;
    }

    Ok(())
}
