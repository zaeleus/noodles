use std::io::{self, Write};

use super::MISSING;
use crate::{
    variant::{record::Filters as _, record_buf::Filters},
    Header,
};

pub(super) fn write_filters<W>(writer: &mut W, header: &Header, filters: &Filters) -> io::Result<()>
where
    W: Write,
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
        fn t(
            buf: &mut Vec<u8>,
            header: &Header,
            filters: &Filters,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_filters(buf, header, filters)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let header = Header::default();
        let mut buf = Vec::new();

        t(&mut buf, &header, &Filters::default(), b".")?;
        t(&mut buf, &header, &Filters::pass(), b"PASS")?;

        let filters = [String::from("q10")].into_iter().collect();
        t(&mut buf, &header, &filters, b"q10")?;

        let filters = [String::from("q10"), String::from("s50")]
            .into_iter()
            .collect();
        t(&mut buf, &header, &filters, b"q10;s50")?;

        Ok(())
    }
}
