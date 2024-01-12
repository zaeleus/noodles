use std::{error, fmt, str};

use crate::Header;

/// An error when a raw SAM record reference sequence ID fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
    /// The reference sequence is not in the reference sequence dictionary.
    MissingReferenceSequenceDictionaryEntry(String),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => write!(f, "invalid input"),
            Self::MissingReferenceSequenceDictionaryEntry(name) => {
                write!(
                    f,
                    "missing reference sequence dictionary entry for '{name}'"
                )
            }
        }
    }
}

pub(super) fn parse_reference_sequence_id(
    header: &Header,
    src: &[u8],
) -> Result<usize, ParseError> {
    str::from_utf8(src)
        .map_err(|_| ParseError::Invalid)
        .and_then(|s| {
            header
                .reference_sequences()
                .get_index_of(s)
                .ok_or_else(|| ParseError::MissingReferenceSequenceDictionaryEntry(s.into()))
        })
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use super::*;
    use crate::header::record::value::{map::ReferenceSequence, Map};

    #[test]
    fn test_parse_reference_sequence_id() -> Result<(), Box<dyn std::error::Error>> {
        let header = Header::builder()
            .add_reference_sequence(
                "sq0".parse()?,
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(8)?),
            )
            .add_reference_sequence(
                "sq1".parse()?,
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?),
            )
            .build();

        assert_eq!(parse_reference_sequence_id(&header, b"sq0")?, 0);
        assert_eq!(parse_reference_sequence_id(&header, b"sq1")?, 1);

        assert_eq!(
            parse_reference_sequence_id(&header, b"\xf0"),
            Err(ParseError::Invalid),
        );

        assert_eq!(
            parse_reference_sequence_id(&header, b"*"),
            Err(ParseError::MissingReferenceSequenceDictionaryEntry(
                String::from("*")
            ))
        );

        assert_eq!(
            parse_reference_sequence_id(&header, b"sq2"),
            Err(ParseError::MissingReferenceSequenceDictionaryEntry(
                String::from("sq2")
            ))
        );

        Ok(())
    }
}
