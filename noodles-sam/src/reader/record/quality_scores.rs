use std::{error, fmt, mem};

use crate::record::{quality_scores, QualityScores};

/// An error when raw SAM record quality scores fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid(quality_scores::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Empty => None,
            Self::Invalid(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::Invalid(_) => write!(f, "invalid input"),
        }
    }
}

pub(crate) fn parse_quality_scores(
    src: &[u8],
    quality_scores: &mut QualityScores,
) -> Result<(), ParseError> {
    const OFFSET: u8 = b'!';

    if src.is_empty() {
        return Err(ParseError::Empty);
    }

    let raw_quality_scores = Vec::from(mem::take(quality_scores));

    let mut raw_scores: Vec<u8> = raw_quality_scores.into_iter().map(u8::from).collect();
    raw_scores.extend(src.iter().map(|n| n.wrapping_sub(OFFSET)));

    *quality_scores = QualityScores::try_from(raw_scores).map_err(ParseError::Invalid)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        let mut quality_scores = QualityScores::default();

        quality_scores.clear();
        parse_quality_scores(b"NDLS", &mut quality_scores)?;
        let expected = QualityScores::try_from(vec![45, 35, 43, 50])?;
        assert_eq!(quality_scores, expected);

        quality_scores.clear();
        assert_eq!(
            parse_quality_scores(b"", &mut quality_scores),
            Err(ParseError::Empty)
        );

        quality_scores.clear();
        assert!(matches!(
            parse_quality_scores(&[0x07], &mut quality_scores),
            Err(ParseError::Invalid(_))
        ));

        Ok(())
    }
}
