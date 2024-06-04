use std::io::{self, Write};

pub(crate) fn write_string<W>(writer: &mut W, s: &str) -> io::Result<()>
where
    W: Write,
{
    // § 1.4 "Meta-information lines" (2024-04-20): "An _unstructured_ meta-information line
    // consists of [...] a _value_ (which may not be empty and must not start with a '<'
    // character)..."
    fn is_valid(s: &str) -> bool {
        const LESS_THAN_SIGN: char = '<';
        s.chars().next().is_some_and(|c| c != LESS_THAN_SIGN)
    }

    if is_valid(s) {
        writer.write_all(s.as_bytes())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid unstructured header record value",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_string() -> io::Result<()> {
        let mut buf = Vec::new();

        write_string(&mut buf, "ndls")?;
        assert_eq!(buf, b"ndls");

        assert!(matches!(
            write_string(&mut buf, ""),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        assert!(matches!(
            write_string(&mut buf, "<"),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
