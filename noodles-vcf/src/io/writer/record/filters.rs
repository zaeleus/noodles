use std::io::{self, Write};

use super::MISSING;
use crate::{variant::record::Filters, Header};

pub(super) fn write_filters<W, F>(writer: &mut W, header: &Header, filters: F) -> io::Result<()>
where
    W: Write,
    F: Filters,
{
    const DELIMITER: &[u8] = b";";

    if filters.is_empty() {
        writer.write_all(MISSING)?;
    } else {
        for (i, result) in filters.iter(header).enumerate() {
            let id = result?;

            if i > 0 {
                writer.write_all(DELIMITER)?;
            }

            if is_valid(id) {
                writer.write_all(id.as_bytes())?;
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "invalid filter",
                ));
            }
        }
    }

    Ok(())
}

// § 1.6.1.7 "Fixed fields: FILTER" (2023-08-23): "...no whitespace or semicolons permitted..."
fn is_valid(s: &str) -> bool {
    fn is_valid_char(c: char) -> bool {
        !c.is_whitespace() && c != ';'
    }

    s.chars().all(is_valid_char)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record_buf::Filters as FiltersBuf;

    #[test]
    fn test_write_filters() -> io::Result<()> {
        fn t(
            buf: &mut Vec<u8>,
            header: &Header,
            filters: &FiltersBuf,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_filters(buf, header, filters)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let header = Header::default();
        let mut buf = Vec::new();

        t(&mut buf, &header, &FiltersBuf::default(), b".")?;
        t(&mut buf, &header, &FiltersBuf::pass(), b"PASS")?;

        let filters = [String::from("q10")].into_iter().collect();
        t(&mut buf, &header, &filters, b"q10")?;

        let filters = [String::from("q10"), String::from("s50")]
            .into_iter()
            .collect();
        t(&mut buf, &header, &filters, b"q10;s50")?;

        buf.clear();
        let filters = [String::from("q 10")].into_iter().collect();
        assert!(matches!(
            write_filters(&mut buf, &header, &filters),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid("PASS"));
        assert!(is_valid("q10"));

        assert!(!is_valid("q 10"));
        assert!(!is_valid("q10;"));
    }
}
