//! SAM record data field character.

use std::{error, fmt};

/// A SAM record data field character value.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct Character(u8);

/// An error returned when a raw character fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => f.write_str("invalid input"),
        }
    }
}

impl TryFrom<char> for Character {
    type Error = ParseError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        if c.is_ascii_graphic() {
            // SAFETY: `c` is guaranteed to be '!'..='~'.
            Ok(Self(c as u8))
        } else {
            Err(ParseError::Invalid)
        }
    }
}

impl TryFrom<u8> for Character {
    type Error = ParseError;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        if n.is_ascii_graphic() {
            Ok(Self(n))
        } else {
            Err(ParseError::Invalid)
        }
    }
}

impl From<Character> for char {
    fn from(character: Character) -> Self {
        char::from(character.0)
    }
}

impl From<Character> for u8 {
    fn from(character: Character) -> Self {
        character.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_u8_for_character() {
        assert_eq!(Character::try_from(b'n'), Ok(Character(b'n')));
        assert_eq!(Character::try_from(b'\t'), Err(ParseError::Invalid));
    }

    #[test]
    fn test_try_from_char_for_character() {
        assert_eq!(Character::try_from('n'), Ok(Character(b'n')));
        assert_eq!(Character::try_from('\t'), Err(ParseError::Invalid));
    }

    #[test]
    fn test_from_character_for_u8() {
        assert_eq!(u8::from(Character(b'n')), b'n');
    }

    #[test]
    fn test_from_character_for_char() {
        assert_eq!(char::from(Character(b'n')), 'n');
    }
}
