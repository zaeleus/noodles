use std::io::{self, Write};

use super::MISSING;
use crate::variant::record::AlternateBases;

pub(super) fn write_alternate_bases<W, B>(writer: &mut W, alternate_bases: B) -> io::Result<()>
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

            if is_valid(bases) {
                write!(writer, "{bases}")?;
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "invalid alternate bases",
                ));
            }
        }
    }

    Ok(())
}

// § 1.6.1.5 "Fixed fields: ALT" (2023-08-23): "...no whitespace, commas..."
fn is_valid(s: &str) -> bool {
    fn is_valid_char(c: char) -> bool {
        !c.is_whitespace() && c != ','
    }

    s.chars().all(is_valid_char)
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

        buf.clear();
        let alternate_bases = AlternateBasesBuf::from(vec![String::from("C,")]);
        assert!(matches!(
            write_alternate_bases(&mut buf, &alternate_bases),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid("C"));
        assert!(is_valid("GT"));
        assert!(is_valid("*"));
        assert!(is_valid("."));
        assert!(is_valid("<DEL>"));
        assert!(is_valid("<*>"));
        assert!(is_valid("<NON_REF>"));
        assert!(is_valid("C[1:1["));
        assert!(is_valid("]1:0]A"));

        assert!(!is_valid("C,"));
        assert!(!is_valid("G T"));
    }
}
