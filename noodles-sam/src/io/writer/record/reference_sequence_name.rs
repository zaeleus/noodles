use std::io::{self, Write};

pub(super) fn write_reference_sequence_name<W>(
    writer: &mut W,
    reference_sequence_name: Option<&[u8]>,
) -> io::Result<()>
where
    W: Write,
{
    use super::MISSING;

    let buf = reference_sequence_name.unwrap_or(&[MISSING]);
    writer.write_all(buf)
}

pub(super) fn write_mate_reference_sequence_name<W>(
    writer: &mut W,
    reference_sequence_name: Option<&[u8]>,
    mate_reference_sequence_name: Option<&[u8]>,
) -> io::Result<()>
where
    W: Write,
{
    const EQ: u8 = b'=';

    match (reference_sequence_name, mate_reference_sequence_name) {
        (Some(n), Some(m)) if n == m => writer.write_all(&[EQ]),
        _ => write_reference_sequence_name(writer, mate_reference_sequence_name),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_reference_sequence_name() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_reference_sequence_name(&mut buf, None)?;
        assert_eq!(buf, b"*");

        buf.clear();
        write_reference_sequence_name(&mut buf, Some(b"sq0"))?;
        assert_eq!(buf, b"sq0");

        Ok(())
    }

    #[test]
    fn test_write_mate_reference_sequence_name() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_mate_reference_sequence_name(&mut buf, None, None)?;
        assert_eq!(buf, b"*");

        buf.clear();
        write_mate_reference_sequence_name(&mut buf, Some(b"sq0"), None)?;
        assert_eq!(buf, b"*");

        buf.clear();
        write_mate_reference_sequence_name(&mut buf, None, Some(b"sq0"))?;
        assert_eq!(buf, b"sq0");

        buf.clear();
        write_mate_reference_sequence_name(&mut buf, Some(b"sq0"), Some(b"sq0"))?;
        assert_eq!(buf, b"=");

        buf.clear();
        write_mate_reference_sequence_name(&mut buf, Some(b"sq0"), Some(b"sq1"))?;
        assert_eq!(buf, b"sq1");

        Ok(())
    }
}
