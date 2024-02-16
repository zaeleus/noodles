use std::io::{self, Write};

use super::MISSING;
use crate::variant::{record::AlternateBases as _, record_buf::AlternateBases};

pub(super) fn write_alternate_bases<W>(
    writer: &mut W,
    alternate_bases: &AlternateBases,
) -> io::Result<()>
where
    W: Write,
{
    if alternate_bases.is_empty() {
        writer.write_all(MISSING)?;
    } else {
        for (i, allele) in alternate_bases.as_ref().iter().enumerate() {
            if i > 0 {
                write!(writer, ",")?;
            }

            write!(writer, "{allele}")?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_alternate_bases() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        let alternate_bases = AlternateBases::default();
        write_alternate_bases(&mut buf, &alternate_bases)?;
        assert_eq!(buf, b".");

        buf.clear();
        let alternate_bases = AlternateBases::from(vec![String::from("C")]);
        write_alternate_bases(&mut buf, &alternate_bases)?;
        assert_eq!(buf, b"C");

        buf.clear();
        let alternate_bases = AlternateBases::from(vec![String::from("C"), String::from("GT")]);
        write_alternate_bases(&mut buf, &alternate_bases)?;
        assert_eq!(buf, b"C,GT");

        Ok(())
    }
}
