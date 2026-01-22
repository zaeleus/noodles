use std::{error, fmt};

use crate::alignment::record_buf::QualityScores;

/// An error when raw SAM record quality scores fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
    /// The length does not match the sequence length.
    LengthMismatch {
        /// The actual length.
        actual: usize,
        /// The expected length.
        expected: usize,
    },
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::Invalid => write!(f, "invalid input"),
            Self::LengthMismatch { actual, expected } => {
                write!(f, "length mismatch: expected {expected}, got {actual}")
            }
        }
    }
}

pub(super) fn parse_quality_scores(
    src: &[u8],
    sequence_len: usize,
    quality_scores: &mut QualityScores,
) -> Result<(), ParseError> {
    const OFFSET: u8 = b'!';

    if src.is_empty() {
        return Err(ParseError::Empty);
    } else if src.len() != sequence_len {
        return Err(ParseError::LengthMismatch {
            actual: src.len(),
            expected: sequence_len,
        });
    } else if !is_valid(src) {
        return Err(ParseError::Invalid);
    }

    quality_scores.as_mut().extend(src.iter().map(|n| {
        // SAFETY: `n` is guaranteed to be [33, 126].
        n - OFFSET
    }));

    Ok(())
}

fn is_valid(scores: &[u8]) -> bool {
    scores.iter().all(|n| n.is_ascii_graphic())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_quality_scores() -> Result<(), ParseError> {
        let mut quality_scores = QualityScores::default();

        quality_scores.as_mut().clear();
        parse_quality_scores(b"NDLS", 4, &mut quality_scores)?;
        let expected = [45, 35, 43, 50].into_iter().collect();
        assert_eq!(quality_scores, expected);

        quality_scores.as_mut().clear();
        assert_eq!(
            parse_quality_scores(b"", 0, &mut quality_scores),
            Err(ParseError::Empty)
        );

        quality_scores.as_mut().clear();
        assert_eq!(
            parse_quality_scores(b"NDLS", 2, &mut quality_scores),
            Err(ParseError::LengthMismatch {
                actual: 4,
                expected: 2
            })
        );

        quality_scores.as_mut().clear();
        assert_eq!(
            parse_quality_scores(&[0x08], 1, &mut quality_scores),
            Err(ParseError::Invalid)
        );

        quality_scores.as_mut().clear();
        assert_eq!(
            parse_quality_scores(&[0x7f], 1, &mut quality_scores),
            Err(ParseError::Invalid)
        );

        Ok(())
    }
}
