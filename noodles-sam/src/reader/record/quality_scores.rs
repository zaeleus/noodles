use std::io;

use crate::record::QualityScores;

pub(crate) fn parse_quality_scores(src: &[u8]) -> io::Result<QualityScores> {
    const MISSING: &[u8] = b"*";
    const OFFSET: u8 = b'!';

    if src == MISSING {
        return Ok(QualityScores::default());
    }

    let mut scores = Vec::new();
    scores.extend(src.iter().map(|n| n.wrapping_sub(OFFSET)));

    QualityScores::try_from(scores).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        let actual = parse_quality_scores(b"")?;
        let expected = QualityScores::default();
        assert_eq!(actual, expected);

        let actual = parse_quality_scores(b"NDLS")?;
        let expected = QualityScores::try_from(vec![45, 35, 43, 50])?;
        assert_eq!(actual, expected);

        assert!(matches!(
            parse_quality_scores(&[0x07]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
