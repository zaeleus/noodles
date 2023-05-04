//! SAM record quality scores score.

use std::{
    error,
    fmt::{self, Write},
};

const START_CHAR: char = '!';
const END_CHAR: char = '~';

const OFFSET: u8 = b'!';

/// A SAM record quality scores score.
///
/// A quality score ranges from 0 to 93 (inclusive), where higher is better.
///
/// Quality scores can be represented as ASCII characters. Each score is offset by 33 (`!`) to only
/// use the set of printable characters (`!`-`~`, excluding the space character).
#[derive(Clone, Copy, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub struct Score(pub(crate) u8);

impl Score {
    /// The minimum score (0).
    pub const MIN: Self = Self(0);

    /// The maximum score (93).
    pub const MAX: Self = Self(93);

    /// Returns the inner value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::quality_scores::Score;
    /// assert_eq!(Score::MIN.get(), 0);
    /// ```
    pub const fn get(&self) -> u8 {
        self.0
    }
}

impl fmt::Display for Score {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_char(char::from(*self))
    }
}

/// An error returned when parsing a raw SAM record quality score fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid(u32),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(actual) => {
                write!(f, "expected <= {}, got {actual}", Score::MAX.get())
            }
        }
    }
}

impl TryFrom<char> for Score {
    type Error = ParseError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            START_CHAR..=END_CHAR => Ok(Self((c as u8) - OFFSET)),
            _ => Err(ParseError::Invalid(u32::from(c))),
        }
    }
}

impl TryFrom<u8> for Score {
    type Error = ParseError;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        if n <= Self::MAX.get() {
            Ok(Self(n))
        } else {
            Err(ParseError::Invalid(u32::from(n)))
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
        Self::from(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_char_for_score() {
        assert_eq!(Score::try_from('N'), Ok(Score(45)));
        assert_eq!(
            Score::try_from(' '),
            Err(ParseError::Invalid(u32::from(' ')))
        );
    }

    #[test]
    fn test_try_from_u8_for_score() {
        assert_eq!(Score::try_from(8), Ok(Score(8)));
        assert_eq!(Score::try_from(144), Err(ParseError::Invalid(144)));
    }

    #[test]
    fn test_from_score_for_u8() {
        assert_eq!(u8::from(Score(8)), 8);
    }

    #[test]
    fn test_from_score_for_char() {
        assert_eq!(char::from(Score(45)), 'N');
    }
}
