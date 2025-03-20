use std::io::{self, Write};

pub(super) fn write_type<W>(writer: &mut W, ty: &str) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(ty.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_type() -> io::Result<()> {
        let mut buf = Vec::new();
        write_type(&mut buf, "exon")?;
        assert_eq!(buf, b"exon");
        Ok(())
    }
}
