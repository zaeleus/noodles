use std::io::{self, Write};

use super::MISSING;
use crate::alignment::record::QualityScores;

pub(super) fn write_quality_scores<W, S>(
    writer: &mut W,
    base_count: usize,
    quality_scores: S,
) -> io::Result<()>
where
    W: Write,
    S: QualityScores,
{
    const OFFSET: u8 = b'!';

    if quality_scores.is_empty() {
        writer.write_all(&[MISSING])?;
    } else if quality_scores.len() == base_count {
        if !is_valid(quality_scores.iter()) {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }

        for score in quality_scores.iter() {
            let n = score + OFFSET;
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

fn is_valid<I>(mut scores: I) -> bool
where
    I: Iterator<Item = u8>,
{
    const MAX_SCORE: u8 = b'~';
    scores.all(|score| score <= MAX_SCORE)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::record_buf::QualityScores as QualityScoresBuf;

    #[test]
    fn test_write_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            buf: &mut Vec<u8>,
            base_count: usize,
            quality_scores: &QualityScoresBuf,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_quality_scores(buf, base_count, quality_scores)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 0, &QualityScoresBuf::default(), &[b'*'])?;
        t(&mut buf, 4, &QualityScoresBuf::default(), &[b'*'])?;

        let quality_scores = QualityScoresBuf::from(vec![45, 35, 43, 50]);
        t(&mut buf, 4, &quality_scores, &[b'N', b'D', b'L', b'S'])?;

        let quality_scores = QualityScoresBuf::from(vec![45, 35, 43, 50]);
        buf.clear();
        assert!(matches!(
            write_quality_scores(&mut buf, 3, &quality_scores),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        let quality_scores = QualityScoresBuf::from(vec![255]);
        buf.clear();
        assert!(matches!(
            write_quality_scores(&mut buf, 1, &quality_scores),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
