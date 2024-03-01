use std::io::{self, Write};

use super::MISSING;
use crate::variant::record::AlternateBases;

pub(super) fn write_alternate_bases<W, B>(writer: &mut W, alternate_bases: &B) -> io::Result<()>
where
    W: Write,
    B: AlternateBases,
{
    if alternate_bases.is_empty() {
        writer.write_all(MISSING)?;
    } else {
        for (i, result) in alternate_bases.iter().enumerate() {
            let bases = result?;

            if i > 0 {
                write!(writer, ",")?;
            }

            write!(writer, "{bases}")?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record_buf::AlternateBases as AlternateBasesBuf;

    #[test]
    fn test_write_alternate_bases() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        let alternate_bases = AlternateBasesBuf::default();
        write_alternate_bases(&mut buf, &alternate_bases)?;
        assert_eq!(buf, b".");

        buf.clear();
        let alternate_bases = AlternateBasesBuf::from(vec![String::from("C")]);
        write_alternate_bases(&mut buf, &alternate_bases)?;
        assert_eq!(buf, b"C");

        buf.clear();
        let alternate_bases = AlternateBasesBuf::from(vec![String::from("C"), String::from("GT")]);
        write_alternate_bases(&mut buf, &alternate_bases)?;
        assert_eq!(buf, b"C,GT");

        Ok(())
    }
}
