use std::io::{self, Write};

use super::MISSING;
use crate::record::Ids;

pub(super) fn write_ids<W>(writer: &mut W, ids: &Ids) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b";";

    if ids.is_empty() {
        writer.write_all(MISSING)?;
    } else {
        for (i, id) in ids.iter().enumerate() {
            if i > 0 {
                writer.write_all(DELIMITER)?;
            }

            writer.write_all(id.as_bytes())?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_ids() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, ids: &Ids, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_ids(buf, ids)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let ids = Ids::default();
        t(&mut buf, &ids, b".")?;

        let ids = "id0".parse()?;
        t(&mut buf, &ids, b"id0")?;

        let ids = "id0;id1".parse()?;
        t(&mut buf, &ids, b"id0;id1")?;

        Ok(())
    }
}
