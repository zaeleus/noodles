use std::io::{self, Write};

use crate::feature::record_buf::Strand;

pub(super) fn write_strand<W>(writer: &mut W, strand: Option<Strand>) -> io::Result<()>
where
    W: Write,
{
    let c = match strand {
        Some(Strand::Forward) => b"+",
        Some(Strand::Reverse) => b"-",
        None => b".",
    };

    writer.write_all(c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_strand() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_strand(&mut buf, None)?;
        assert_eq!(buf, b".");

        buf.clear();
        write_strand(&mut buf, Some(Strand::Forward))?;
        assert_eq!(buf, b"+");

        buf.clear();
        write_strand(&mut buf, Some(Strand::Reverse))?;
        assert_eq!(buf, b"-");

        Ok(())
    }
}
