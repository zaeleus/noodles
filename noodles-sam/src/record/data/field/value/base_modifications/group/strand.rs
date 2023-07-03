use std::{error, fmt};

/// The strand on which the modification was observed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Strand {
    /// The same strand.
    Forward,
    /// The opposite strand.
    Reverse,
}

/// An error returned when a base modifications group strand fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => write!(f, "invalid input"),
        }
    }
}

impl TryFrom<u8> for Strand {
    type Error = ParseError;

    fn try_from(b: u8) -> Result<Self, Self::Error> {
        match b {
            b'+' => Ok(Self::Forward),
            b'-' => Ok(Self::Reverse),
            _ => Err(ParseError::Invalid),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_u8_for_strand() {
        assert_eq!(Strand::try_from(b'+'), Ok(Strand::Forward));
        assert_eq!(Strand::try_from(b'-'), Ok(Strand::Reverse));

        assert_eq!(Strand::try_from(b'n'), Err(ParseError::Invalid));
    }
}
