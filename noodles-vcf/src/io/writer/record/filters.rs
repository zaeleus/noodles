use std::io::{self, Write};

use super::MISSING;
use crate::variant::record_buf::Filters;

pub(super) fn write_filters<W>(writer: &mut W, filters: Option<&Filters>) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b";";
    const PASS: &[u8] = b"PASS";

    if let Some(filters) = filters {
        match filters {
            Filters::Pass => writer.write_all(PASS)?,
            Filters::Fail(ids) => {
                for (i, id) in ids.iter().enumerate() {
                    if i > 0 {
                        writer.write_all(DELIMITER)?;
                    }

                    writer.write_all(id.as_bytes())?;
                }
            }
        }
    } else {
        writer.write_all(MISSING)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_filters() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, filters: Option<&Filters>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_filters(buf, filters)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, b".")?;
        t(&mut buf, Some(&Filters::Pass), b"PASS")?;

        let filters = Filters::try_from_iter(["q10"])?;
        t(&mut buf, Some(&filters), b"q10")?;

        let filters = Filters::try_from_iter(["q10", "s50"])?;
        t(&mut buf, Some(&filters), b"q10;s50")?;

        Ok(())
    }
}
