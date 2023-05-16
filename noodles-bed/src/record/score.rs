//! BED record score.

use std::{error, fmt, num, str::FromStr};

const MIN_VALUE: u16 = 1;
const MAX_VALUE: u16 = 1000;

/// A BED record score.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Score(u16);

impl Score {
    /// The minimun score (1).
    pub const MIN: Self = Self(MIN_VALUE);

    /// The maximum score (1000).
    pub const MAX: Self = Self(MAX_VALUE);

    /// Creates a score if the given value is within range.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::record::Score;
    ///
    /// assert!(Score::new(8).is_some());
    ///
    /// assert!(Score::new(0).is_none());
    /// assert!(Score::new(1001).is_none());
    /// ```
    pub const fn new(n: u16) -> Option<Self> {
        if MIN_VALUE <= n && n <= MAX_VALUE {
            Some(Self(n))
        } else {
            None
        }
    }

    /// Returns the inner value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::record::Score;
    /// let score = Score::new(8).unwrap();
    /// assert_eq!(score.get(), 8);
    /// ```
    pub const fn get(&self) -> u16 {
        self.0
    }
}

impl fmt::Display for Score {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// An error returned when a raw BED record score fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input failed to be parsed as an integer.
    Parse(num::ParseIntError),
    /// The input is invalid.
    Invalid(TryFromIntError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Parse(e) => Some(e),
            Self::Invalid(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(_) => f.write_str("parse error"),
            Self::Invalid(_) => f.write_str("invalid input"),
        }
    }
}

impl FromStr for Score {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let n: u16 = s.parse().map_err(ParseError::Parse)?;
        Self::try_from(n).map_err(ParseError::Invalid)
    }
}

/// An error returned when a raw BED record score fails to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromIntError(u16);

impl error::Error for TryFromIntError {}

impl fmt::Display for TryFromIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid value: {}", self.0)
    }
}

impl TryFrom<u16> for Score {
    type Error = TryFromIntError;

    fn try_from(n: u16) -> Result<Self, Self::Error> {
        if (MIN_VALUE..=MAX_VALUE).contains(&n) {
            Ok(Self(n))
        } else {
            Err(TryFromIntError(n))
        }
    }
}

impl From<Score> for u16 {
    fn from(score: Score) -> Self {
        score.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Score(1).to_string(), "1");
        assert_eq!(Score(1000).to_string(), "1000");
    }

    #[test]
    fn test_try_from_u16_for_score() {
        assert_eq!(Score::try_from(1), Ok(Score::MIN));
        assert_eq!(Score::try_from(8), Ok(Score(8)));
        assert_eq!(Score::try_from(1000), Ok(Score::MAX));

        assert_eq!(Score::try_from(0), Err(TryFromIntError(0)));
        assert_eq!(Score::try_from(1001), Err(TryFromIntError(1001)));
    }

    #[test]
    fn test_from_score_for_u16() {
        assert_eq!(u16::from(Score(8)), 8);
    }
}
