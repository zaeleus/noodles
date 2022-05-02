use std::io::{self, Write};

use super::MISSING;
use crate::alignment::record::AlignmentQualityScores;

pub fn write_quality_scores<W, S>(writer: &mut W, quality_scores: &S) -> io::Result<()>
where
    W: Write,
    S: AlignmentQualityScores + ?Sized,
{
    const MIN_VALUE: u8 = b'!';

    if quality_scores.is_empty() {
        writer.write_all(&[MISSING])?;
    } else {
        for score in quality_scores.scores() {
            let n = u8::from(score) + MIN_VALUE;
            writer.write_all(&[n])?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        use crate::alignment::record::QualityScores;

        let mut buf = Vec::new();
        write_quality_scores(&mut buf, &QualityScores::default())?;
        assert_eq!(buf, b"*");

        buf.clear();
        let quality_scores: QualityScores = "NDLS".parse()?;
        write_quality_scores(&mut buf, &quality_scores)?;
        assert_eq!(buf, b"NDLS");

        Ok(())
    }
}
