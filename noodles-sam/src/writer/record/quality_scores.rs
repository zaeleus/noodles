use std::io::{self, Write};

use super::MISSING;
use crate::record::QualityScores;

pub fn write_quality_scores<W>(
    writer: &mut W,
    base_count: usize,
    quality_scores: &QualityScores,
) -> io::Result<()>
where
    W: Write,
{
    const MIN_VALUE: u8 = b'!';

    if quality_scores.is_empty() {
        writer.write_all(&[MISSING])?;
    } else if quality_scores.len() == base_count {
        for &score in quality_scores.as_ref() {
            let n = u8::from(score) + MIN_VALUE;
            writer.write_all(&[n])?;
        }
    } else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "sequence-quality scores length mismatch: expected {}, got {}",
                base_count,
                quality_scores.len()
            ),
        ));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            buf: &mut Vec<u8>,
            base_count: usize,
            quality_scores: &QualityScores,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_quality_scores(buf, base_count, quality_scores)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 0, &QualityScores::default(), &[b'*'])?;
        t(&mut buf, 4, &QualityScores::default(), &[b'*'])?;
        t(&mut buf, 4, &"NDLS".parse()?, &[b'N', b'D', b'L', b'S'])?;

        let quality_scores = "NDLS".parse()?;
        buf.clear();
        assert!(matches!(
            write_quality_scores(&mut buf, 3, &quality_scores),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
