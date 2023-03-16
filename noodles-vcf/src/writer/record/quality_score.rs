use std::io::{self, Write};

use super::MISSING;
use crate::record::QualityScore;

pub(super) fn write_quality_score<W>(
    writer: &mut W,
    quality_score: Option<QualityScore>,
) -> io::Result<()>
where
    W: Write,
{
    if let Some(qual) = quality_score.map(f32::from) {
        write!(writer, "{qual}")?;
    } else {
        writer.write_all(MISSING)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_quality_score() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            buf: &mut Vec<u8>,
            quality_score: Option<QualityScore>,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_quality_score(buf, quality_score)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, b".")?;
        t(&mut buf, QualityScore::try_from(0.0).map(Some)?, b"0")?;
        t(&mut buf, QualityScore::try_from(8.13).map(Some)?, b"8.13")?;

        Ok(())
    }
}
