use std::io::{self, Write};

use crate::variant::record::samples::keys::key;

pub(super) fn write_keys<'a, W, I>(writer: &mut W, keys: I) -> io::Result<()>
where
    W: Write,
    I: Iterator<Item = io::Result<&'a str>>,
{
    const DELIMITER: &[u8] = b":";

    for (i, result) in keys.enumerate() {
        let key = result?;

        if i > 0 {
            if key == key::GENOTYPE {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "GT must be first series",
                ));
            }

            writer.write_all(DELIMITER)?;
        }

        writer.write_all(key.as_bytes())?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_keys() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        let keys = [Ok("GT")];
        write_keys(&mut buf, keys.into_iter())?;
        assert_eq!(buf, b"GT");

        buf.clear();
        let keys = [Ok("GT"), Ok("GQ")];
        write_keys(&mut buf, keys.into_iter())?;
        assert_eq!(buf, b"GT:GQ");

        buf.clear();
        let keys = [Ok("GQ")];
        write_keys(&mut buf, keys.into_iter())?;
        assert_eq!(buf, b"GQ");

        buf.clear();
        let keys = [Ok("GQ"), Ok("GT")];
        assert!(matches!(
            write_keys(&mut buf, keys.into_iter()),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
