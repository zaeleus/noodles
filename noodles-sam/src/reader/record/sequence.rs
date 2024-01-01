use std::{error, fmt};

use crate::alignment::record_buf::Sequence;

/// An error when a raw SAM record sequence fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
        }
    }
}

pub(super) fn parse_sequence(src: &[u8], sequence: &mut Sequence) -> Result<(), ParseError> {
    if src.is_empty() {
        return Err(ParseError::Empty);
    }

    let sequence = sequence.as_mut();
    sequence.extend(src);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sequence() -> Result<(), ParseError> {
        let mut sequence = Sequence::default();

        sequence.as_mut().clear();
        parse_sequence(b"ACGT", &mut sequence)?;
        let expected = Sequence::from(b"ACGT".to_vec());
        assert_eq!(sequence, expected);

        sequence.as_mut().clear();
        parse_sequence(&[0x07], &mut sequence)?;
        let expected = Sequence::from(vec![0x07]);
        assert_eq!(sequence, expected);

        sequence.as_mut().clear();
        assert_eq!(parse_sequence(b"", &mut sequence), Err(ParseError::Empty));

        Ok(())
    }
}
