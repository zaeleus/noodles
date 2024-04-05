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

            if is_valid(id) {
                writer.write_all(id.as_bytes())?;
            } else {
                return Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid ID"));
            }
        }
    }

    Ok(())
}

// § 1.6.1.3 "Fixed fields: ID" (2023-08-23): "...no whitespace or semicolons permitted..."
fn is_valid(s: &str) -> bool {
    fn is_valid_char(c: char) -> bool {
        !c.is_whitespace() && c != ';'
    }

    s.chars().all(is_valid_char)
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

        buf.clear();
        let ids = [String::from("id 0")].into_iter().collect();
        assert!(matches!(
            write_ids(&mut buf, &ids),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid("id0"));

        assert!(!is_valid("id 0"));
        assert!(!is_valid("id0;"));
    }
}
