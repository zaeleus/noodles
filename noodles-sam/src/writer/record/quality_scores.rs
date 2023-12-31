use std::io::{self, Write};

use super::MISSING;
use crate::alignment::record_buf::QualityScores;

pub fn write_quality_scores<W>(
    writer: &mut W,
    base_count: usize,
    quality_scores: &QualityScores,
) -> io::Result<()>
where
    W: Write,
{
    const OFFSET: u8 = b'!';

    let quality_scores = quality_scores.as_ref();

    if quality_scores.is_empty() {
        writer.write_all(&[MISSING])?;
    } else if quality_scores.len() == base_count {
        if !is_valid(quality_scores) {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }

        for &score in quality_scores {
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

fn is_valid(scores: &[u8]) -> bool {
    const MAX_SCORE: u8 = b'~';
    scores.iter().all(|&score| score <= MAX_SCORE)
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

        let quality_scores = QualityScores::from(vec![45, 35, 43, 50]);
        t(&mut buf, 4, &quality_scores, &[b'N', b'D', b'L', b'S'])?;

        let quality_scores = QualityScores::from(vec![45, 35, 43, 50]);
        buf.clear();
        assert!(matches!(
            write_quality_scores(&mut buf, 3, &quality_scores),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        let quality_scores = QualityScores::from(vec![255]);
        buf.clear();
        assert!(matches!(
            write_quality_scores(&mut buf, 1, &quality_scores),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
