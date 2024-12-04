use std::io::{self, Write};

pub(super) fn write_reference_sequence_name<W>(
    writer: &mut W,
    reference_sequence_name: &str,
) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(reference_sequence_name.as_bytes())
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

        Ok(())
    }
}
