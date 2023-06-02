use std::{error, fmt, mem};

use bytes::Buf;

/// An error when a raw BAM record reference sequence ID fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The input is invalid.
    Invalid,
    /// The reference sequence is not in the reference sequence dictionary.
    MissingReferenceSequenceDictionaryEntry { actual: usize, expected: usize },
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::Invalid => write!(f, "invalid input"),
            Self::MissingReferenceSequenceDictionaryEntry { actual, expected } => {
                write!(
                    f,
                    "missing reference sequence dictionary entry: expected id < {expected}, got {actual}"
                )
            }
        }
    }
}

pub(crate) fn get_reference_sequence_id<B>(
    src: &mut B,
    n_ref: usize,
) -> Result<Option<usize>, ParseError>
where
    B: Buf,
{
    const UNMAPPED: i32 = -1;

    if src.remaining() < mem::size_of::<i32>() {
        return Err(ParseError::UnexpectedEof);
    }

    match src.get_i32_le() {
        UNMAPPED => Ok(None),
        n => usize::try_from(n)
            .map_err(|_| ParseError::Invalid)
            .and_then(|m| {
                if m < n_ref {
                    Ok(Some(m))
                } else {
                    Err(ParseError::MissingReferenceSequenceDictionaryEntry {
                        actual: m,
                        expected: n_ref,
                    })
                }
            }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_reference_sequence_id() {
        let data = (-1i32).to_le_bytes();
        let mut src = &data[..];
        assert_eq!(get_reference_sequence_id(&mut src, 1), Ok(None));

        let data = 0i32.to_le_bytes();
        let mut src = &data[..];
        assert_eq!(get_reference_sequence_id(&mut src, 1), Ok(Some(0)));

        let data = [];
        let mut src = &data[..];
        assert_eq!(
            get_reference_sequence_id(&mut src, 1),
            Err(ParseError::UnexpectedEof)
        );

        let data = (-2i32).to_le_bytes();
        let mut src = &data[..];
        assert_eq!(
            get_reference_sequence_id(&mut src, 1),
            Err(ParseError::Invalid)
        );

        let data = 0i32.to_le_bytes();
        let mut src = &data[..];
        assert_eq!(
            get_reference_sequence_id(&mut src, 0),
            Err(ParseError::MissingReferenceSequenceDictionaryEntry {
                actual: 0,
                expected: 0
            })
        );
    }
}
