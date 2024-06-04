use std::io::{self, Write};

use crate::header::FileFormat;

const VCF_4_3: FileFormat = FileFormat::new(4, 3);

pub(crate) fn write_string<W>(writer: &mut W, file_format: FileFormat, s: &str) -> io::Result<()>
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

    if file_format < VCF_4_3 || is_valid(s) {
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
        const VCF_4_2: FileFormat = FileFormat::new(4, 2);

        let mut buf = Vec::new();

        buf.clear();
        write_string(&mut buf, VCF_4_3, "ndls")?;
        assert_eq!(buf, b"ndls");

        buf.clear();
        assert!(matches!(
            write_string(&mut buf, VCF_4_3, ""),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        buf.clear();
        write_string(&mut buf, VCF_4_2, "<")?;
        assert_eq!(buf, b"<");

        buf.clear();
        assert!(matches!(
            write_string(&mut buf, VCF_4_3, "<"),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
