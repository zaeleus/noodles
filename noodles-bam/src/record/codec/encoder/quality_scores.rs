use noodles_sam::alignment::record::QualityScores;
use std::{io, iter};

use super::num::write_u8;

pub fn write_quality_scores<S>(
    dst: &mut Vec<u8>,
    base_count: usize,
    quality_scores: S,
) -> io::Result<()>
where
    S: QualityScores,
{
    // § 4.2.3 SEQ and QUAL encoding (2022-08-22)
    const MISSING: u8 = 255;

    if quality_scores.len() == base_count {
        for result in quality_scores.iter() {
            let n = result?;

            if is_valid_score(n) {
                write_u8(dst, n);
            } else {
                return Err(io::Error::from(io::ErrorKind::InvalidInput));
            }
        }
    } else if quality_scores.is_empty() {
        dst.extend(iter::repeat_n(MISSING, base_count));
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

fn is_valid_score(score: u8) -> bool {
    // § 4.2.3 "SEQ and QUAL encoding" (2023-05-24): "Base qualities are stored as bytes in the
    // range [0, 93]..."
    const MAX_SCORE: u8 = 93;
    score <= MAX_SCORE
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record_buf::QualityScores as QualityScoresBuf;

    use super::*;

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

        t(&mut buf, 0, &QualityScoresBuf::default(), &[])?;
        t(
            &mut buf,
            4,
            &QualityScoresBuf::default(),
            &[0xff, 0xff, 0xff, 0xff],
        )?;

        let quality_scores = QualityScoresBuf::from(vec![45, 35, 43, 50]);
        t(&mut buf, 4, &quality_scores, &[0x2d, 0x23, 0x2b, 0x32])?;

        let quality_scores = QualityScoresBuf::from(vec![45, 35, 43, 50]);
        buf.clear();
        assert!(matches!(
            write_quality_scores(&mut buf, 3, &quality_scores),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
