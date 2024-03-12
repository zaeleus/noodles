use std::io::{self, Write};

use crate::variant::record::ReferenceBases;

pub(super) fn write_reference_bases<W, B>(writer: &mut W, reference_bases: B) -> io::Result<()>
where
    W: Write,
    B: ReferenceBases,
{
    for result in reference_bases.iter() {
        let base = result?;
        writer.write_all(&[base])?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_reference_bases() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_reference_bases(&mut buf, "A")?;
        assert_eq!(buf, b"A");

        buf.clear();
        write_reference_bases(&mut buf, "AC")?;
        assert_eq!(buf, b"AC");

        Ok(())
    }
}
