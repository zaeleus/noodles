use std::io::{self, Write};

use super::MISSING;
use crate::variant::{record::Filters as _, record_buf::Filters};

pub(super) fn write_filters<W>(writer: &mut W, filters: &Filters) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b";";

    if filters.is_empty() {
        writer.write_all(MISSING)?;
    } else {
        for (i, id) in filters.iter().enumerate() {
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
    fn test_write_filters() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, filters: &Filters, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_filters(buf, filters)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Filters::default(), b".")?;
        t(&mut buf, &Filters::pass(), b"PASS")?;

        let filters = [String::from("q10")].into_iter().collect();
        t(&mut buf, &filters, b"q10")?;

        let filters = [String::from("q10"), String::from("s50")]
            .into_iter()
            .collect();
        t(&mut buf, &filters, b"q10;s50")?;

        Ok(())
    }
}
