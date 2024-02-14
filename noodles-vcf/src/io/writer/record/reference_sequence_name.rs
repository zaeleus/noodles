use std::io::{self, Write};

pub(super) fn write_reference_sequence_name<W>(
    writer: &mut W,
    reference_sequence_name: &str,
) -> io::Result<()>
where
    W: Write,
{
    const SYMBOL_PREFIX: char = '<';
    const SYMBOL_SUFFIX: char = '>';

    if let Some(s) = reference_sequence_name.strip_prefix(SYMBOL_PREFIX) {
        if let Some(s) = s.strip_suffix(SYMBOL_SUFFIX) {
            if is_valid(s) {
                return writer.write_all(reference_sequence_name.as_bytes());
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "invalid symbol",
                ));
            }
        }
    }

    if is_valid(reference_sequence_name) {
        writer.write_all(reference_sequence_name.as_bytes())
    } else {
        Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid name"))
    }
}

// § 1.4.7 "Contig field format" (2023-08-23): "`[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*`".
fn is_valid(s: &str) -> bool {
    fn is_valid_char(c: char) -> bool {
        ('!'..='~').contains(&c)
            && !matches!(
                c,
                '\\' | ',' | '"' | '`' | '\'' | '(' | ')' | '[' | ']' | '{' | '}' | '<' | '>',
            )
    }

    let mut chars = s.chars();

    let is_valid_first_char = chars
        .next()
        .map(|c| c != '*' && c != '=' && is_valid_char(c))
        .unwrap_or_default();

    is_valid_first_char && chars.all(is_valid_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_reference_sequence_name() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, reference_sequence_name: &str, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_reference_sequence_name(buf, reference_sequence_name)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, "sq0", b"sq0")?;
        t(&mut buf, "<sq0>", b"<sq0>")?;

        buf.clear();
        assert!(matches!(
            write_reference_sequence_name(&mut buf, "sq 0"),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        buf.clear();
        assert!(matches!(
            write_reference_sequence_name(&mut buf, "<sq 0>"),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid("sq0"));

        assert!(!is_valid(""));
        assert!(!is_valid("sq 0"));
        assert!(!is_valid("sq[0]"));
        assert!(!is_valid(">sq0"));
        assert!(!is_valid("*sq0"));
        assert!(!is_valid("=sq0"));
    }
}
