use std::{io, mem};

use crate::record::QualityScores;

pub(crate) fn parse_quality_scores(
    src: &[u8],
    quality_scores: &mut QualityScores,
) -> io::Result<()> {
    const OFFSET: u8 = b'!';

    let raw_quality_scores = Vec::from(mem::take(quality_scores));

    let mut raw_scores: Vec<u8> = raw_quality_scores.into_iter().map(u8::from).collect();
    raw_scores.extend(src.iter().map(|n| n.wrapping_sub(OFFSET)));

    *quality_scores = QualityScores::try_from(raw_scores)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        let mut quality_scores = QualityScores::default();

        quality_scores.clear();
        parse_quality_scores(b"", &mut quality_scores)?;
        assert!(quality_scores.is_empty());

        quality_scores.clear();
        parse_quality_scores(b"NDLS", &mut quality_scores)?;
        let expected = QualityScores::try_from(vec![45, 35, 43, 50])?;
        assert_eq!(quality_scores, expected);

        quality_scores.clear();
        assert!(matches!(
            parse_quality_scores(&[0x07], &mut quality_scores),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
