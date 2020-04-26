use std::{convert::TryFrom, error, fmt};

const OFFSET: u8 = b'!';

const START_CHAR: char = '!';
const END_CHAR: char = '~';

const MAX_SCORE: u8 = 93;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Score(u8);

#[derive(Debug)]
pub struct TryFromCharError(char);

impl error::Error for TryFromCharError {}

impl fmt::Display for TryFromCharError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
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

#[derive(Debug)]
pub struct TryFromUByteError(u8);

impl error::Error for TryFromUByteError {}

impl fmt::Display for TryFromUByteError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
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
