use std::io::{self, Write};

pub(super) fn write_keys<'a, W, I>(writer: &mut W, keys: I) -> io::Result<()>
where
    W: Write,
    I: Iterator<Item = io::Result<&'a str>>,
{
    const DELIMITER: &[u8] = b":";

    for (i, result) in keys.enumerate() {
        let key = result?;

        if i > 0 {
            writer.write_all(DELIMITER)?;
        }

        writer.write_all(key.as_bytes())?;
    }

    Ok(())
}
