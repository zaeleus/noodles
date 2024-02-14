use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use crate::record::codec::value::Float;

pub(super) fn write_quality_score<W>(writer: &mut W, quality_score: Option<f32>) -> io::Result<()>
where
    W: Write,
{
    let float = quality_score.map(Float::from).unwrap_or(Float::Missing);
    writer.write_f32::<LittleEndian>(f32::from(float))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_qual() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, quality_score: Option<f32>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_quality_score(buf, quality_score)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, &[0x01, 0x00, 0x80, 0x7f])?;
        t(&mut buf, Some(0.0), &[0x00, 0x00, 0x00, 0x00])?;
        t(&mut buf, Some(8.0), &[0x00, 0x00, 0x00, 0x41])?;

        Ok(())
    }
}
