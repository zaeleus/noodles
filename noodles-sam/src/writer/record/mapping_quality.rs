use std::io::{self, Write};

use crate::{alignment::record_buf::MappingQuality, writer::num};

pub(super) fn write_mapping_quality<W>(
    writer: &mut W,
    mapping_quality: Option<MappingQuality>,
) -> io::Result<()>
where
    W: Write,
{
    const MISSING: u8 = 255;

    let n = mapping_quality.map(u8::from).unwrap_or(MISSING);
    num::write_u8(writer, n)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_mapping_quality() -> io::Result<()> {
        fn t(
            buf: &mut Vec<u8>,
            mapping_quality: Option<MappingQuality>,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_mapping_quality(buf, mapping_quality)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, b"255")?;
        t(&mut buf, MappingQuality::new(8), b"8")?;

        Ok(())
    }
}
