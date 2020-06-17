use std::{convert::TryFrom, error, fmt};

const OFFSET: u8 = b'!';

const START_CHAR: char = '!';
const END_CHAR: char = '~';

const MAX_SCORE: u8 = 93;

/// A SAM record quality scores score.
///
/// A quality score ranges from 0 to 93 (inclusive), where higher is better.
///
/// Quality scores can be represented as ASCII characters. Each score is offset by 33 (`!`) to only
/// use the set of printable characters (`!`-`~`, excluding the space character).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Score(u8);

impl fmt::Display for Score {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", char::from(*self))
    }
}

/// An error returned when the conversion from a character to a SAM quality scores score fails.
#[derive(Debug)]
pub struct TryFromCharError(char);

impl error::Error for TryFromCharError {}

impl fmt::Display for TryFromCharError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid score: expected {{{}..={}}}, got {}",
            START_CHAR, END_CHAR, self.0
        )
    }
}

impl TryFrom<char> for Score {
    type Error = TryFromCharError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            START_CHAR..=END_CHAR => Ok(Score((c as u8) - OFFSET)),
            _ => Err(TryFromCharError(c)),
        }
    }
}

/// An error returned when the conversion from a byte to a SAM quality scores score fails.
#[derive(Debug)]
pub struct TryFromUByteError(u8);

impl error::Error for TryFromUByteError {}

impl fmt::Display for TryFromUByteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid score: expected {{0..={}}}, got {}",
            MAX_SCORE, self.0
        )
    }
}

impl TryFrom<u8> for Score {
    type Error = TryFromUByteError;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        if n <= MAX_SCORE {
            Ok(Self(n))
        } else {
            Err(TryFromUByteError(n))
        }
    }
}

impl From<Score> for u8 {
    fn from(score: Score) -> Self {
        score.0
    }
}

impl From<Score> for char {
    fn from(score: Score) -> Self {
        let value = u8::from(score) + OFFSET;
        char::from(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_char_for_score() -> Result<(), TryFromCharError> {
        assert_eq!(Score::try_from('N').map(u8::from)?, 45);
        assert!(Score::try_from(' ').is_err());
        Ok(())
    }

    #[test]
    fn test_try_from_u8_for_score() -> Result<(), TryFromUByteError> {
        assert_eq!(Score::try_from(8).map(u8::from)?, 8);
        assert!(Score::try_from(144).is_err());
        Ok(())
    }

    #[test]
    fn test_from_score_for_u8() -> Result<(), TryFromUByteError> {
        assert_eq!(Score::try_from(8).map(u8::from)?, 8);
        Ok(())
    }

    #[test]
    fn test_from_score_for_char() -> Result<(), TryFromCharError> {
        assert_eq!(Score::try_from('N').map(char::from)?, 'N');
        Ok(())
    }
}
