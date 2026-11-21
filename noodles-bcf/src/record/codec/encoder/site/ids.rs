use std::io::{self, Write};

use noodles_vcf::variant::record::Ids;

use crate::record::codec::{Value, encoder::value::write_value};

const SEPARATOR: char = ';';

pub(super) fn write_ids<W, I>(writer: &mut W, ids: I) -> io::Result<()>
where
    W: Write,
    I: Ids,
{
    if ids.is_empty() {
        return write_value(writer, Some(Value::String(None)));
    }

    let mut s = String::new();

    for (i, id) in ids.iter().enumerate() {
        if !is_valid(id) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("invalid ID: {id:?}"),
            ));
        }

        if i > 0 {
            s.push(SEPARATOR);
        }

        s.push_str(id);
    }

    write_value(writer, Some(Value::String(Some(&s))))
}

// § 1.6.1.3 "Fixed fields: ID" (2026-02-25): "...no whitespace or semicolons permitted..."
fn is_valid(s: &str) -> bool {
    const MISSING: &str = ".";

    fn is_valid_char(c: char) -> bool {
        !c.is_whitespace() && !matches!(c, SEPARATOR | '\0')
    }

    !s.is_empty() && s != MISSING && s.chars().all(is_valid_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_ids() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_vcf::variant::record_buf::Ids as IdsBuf;

        fn t(buf: &mut Vec<u8>, ids: &IdsBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_ids(buf, ids)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &IdsBuf::default(), &[0x07])?;
        t(
            &mut buf,
            &[String::from("nd0")].into_iter().collect(),
            &[0x37, b'n', b'd', b'0'],
        )?;
        t(
            &mut buf,
            &[String::from("nd0"), String::from("nd1")]
                .into_iter()
                .collect(),
            &[0x77, b'n', b'd', b'0', b';', b'n', b'd', b'1'],
        )?;

        buf.clear();
        let ids = [String::new()].into_iter().collect();
        assert!(matches!(
            write_ids(&mut buf, &ids),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid("nd0"));

        assert!(!is_valid(""));
        assert!(!is_valid("."));
        assert!(!is_valid("nd 0"));
        assert!(!is_valid("nd0;"));
        assert!(!is_valid("nd\0"));
    }
}
