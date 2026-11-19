use std::io::{self, Write};

use crate::{io::writer::num::write_f32_le, record::codec::value::Float};

pub(super) fn write_quality_score<W>(writer: &mut W, quality_score: Option<f32>) -> io::Result<()>
where
    W: Write,
{
    let n = match quality_score {
        Some(n) => match Float::from(n) {
            Float::Value(m) => m,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid quality score: {n}"),
                ));
            }
        },
        None => f32::from(Float::Missing),
    };

    write_f32_le(writer, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_qual() -> io::Result<()> {
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

        for float in [
            Float::Missing,
            Float::EndOfVector,
            Float::Reserved(f32::from_bits(0x7f800003)),
        ] {
            buf.clear();

            assert!(matches!(
                write_quality_score(&mut buf, Some(f32::from(float))),
                Err(e) if e.kind() == io::ErrorKind::InvalidInput
            ));
        }

        Ok(())
    }
}
