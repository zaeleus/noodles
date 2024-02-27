use std::io::{self, Write};

use super::MISSING;
use crate::variant::record::Ids;

pub(super) fn write_ids<W, I>(writer: &mut W, ids: I) -> io::Result<()>
where
    W: Write,
    I: Ids,
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
    use crate::variant::record_buf::Ids as IdsBuf;

    #[test]
    fn test_write_ids() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, ids: &IdsBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_ids(buf, ids)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let ids = IdsBuf::default();
        t(&mut buf, &ids, b".")?;

        let ids = [String::from("id0")].into_iter().collect();
        t(&mut buf, &ids, b"id0")?;

        let ids = [String::from("id0"), String::from("id1")]
            .into_iter()
            .collect();
        t(&mut buf, &ids, b"id0;id1")?;

        Ok(())
    }
}
