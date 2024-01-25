use std::io::{self, Write};

use crate::Header;

pub(super) fn write_header<W>(writer: &mut W, header: &Header) -> io::Result<()>
where
    W: Write,
{
    write!(writer, "{header}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_header() -> io::Result<()> {
        let mut buf = Vec::new();

        let header = Header::default();
        write_header(&mut buf, &header)?;

        let expected = b"##fileformat=VCFv4.4
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(buf, expected);

        Ok(())
    }
}
