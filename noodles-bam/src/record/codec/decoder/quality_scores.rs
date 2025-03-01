use std::{error, fmt};

use noodles_sam::alignment::record_buf::QualityScores;

/// An error when raw BAM record quality scores fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The input is invalid.
    Invalid,
}

impl error::Error for DecodeError {}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::Invalid => write!(f, "invalid input"),
        }
    }
}

pub fn read_quality_scores(
    src: &mut &[u8],
    quality_scores: &mut QualityScores,
    base_count: usize,
) -> Result<(), DecodeError> {
    let dst = quality_scores.as_mut();

    if base_count == 0 {
        dst.clear();
        return Ok(());
    }

    let (buf, rest) = src
        .split_at_checked(base_count)
        .ok_or(DecodeError::UnexpectedEof)?;

    *src = rest;

    if is_missing_quality_scores(buf) {
        dst.clear();
    } else {
        dst.resize(base_count, 0);
        dst.copy_from_slice(buf);
    }

    Ok(())
}

fn is_missing_quality_scores(src: &[u8]) -> bool {
    const MISSING: u8 = 0xff;

    src.iter().all(|&b| b == MISSING)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_quality_scores() -> Result<(), DecodeError> {
        fn t(mut src: &[u8], expected: &QualityScores) -> Result<(), DecodeError> {
            let mut actual = QualityScores::default();
            read_quality_scores(&mut src, &mut actual, expected.as_ref().len())?;
            assert_eq!(&actual, expected);
            Ok(())
        }

        t(&[], &QualityScores::default())?;
        t(
            &[0x2d, 0x23, 0x2b, 0x32],
            &QualityScores::from(vec![45, 35, 43, 50]),
        )?;

        Ok(())
    }

    #[test]
    fn test_get_quality_scores_with_sequence_and_no_quality_scores() -> Result<(), DecodeError> {
        let data = [0xff, 0xff, 0xff, 0xff];
        let mut buf = &data[..];

        let mut quality_scores = QualityScores::default();
        read_quality_scores(&mut buf, &mut quality_scores, 4)?;

        assert!(quality_scores.is_empty());

        Ok(())
    }
}
