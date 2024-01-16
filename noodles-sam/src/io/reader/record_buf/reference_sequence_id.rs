use std::{error, fmt, str};

use crate::Header;

/// An error when a raw SAM record reference sequence ID fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The reference sequence is not in the reference sequence dictionary.
    MissingReferenceSequenceDictionaryEntry(Vec<u8>),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingReferenceSequenceDictionaryEntry(buf) => {
                let name = str::from_utf8(buf).map_err(|_| fmt::Error)?;

                write!(f, "missing reference sequence dictionary entry: {name}")
            }
        }
    }
}

pub(super) fn parse_reference_sequence_id(
    header: &Header,
    src: &[u8],
) -> Result<usize, ParseError> {
    header
        .reference_sequences()
        .get_index_of(src)
        .ok_or_else(|| ParseError::MissingReferenceSequenceDictionaryEntry(src.into()))
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use super::*;
    use crate::header::record::value::{map::ReferenceSequence, Map};

    #[test]
    fn test_parse_reference_sequence_id() {
        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        const SQ1_LN: NonZeroUsize = match NonZeroUsize::new(13) {
            Some(length) => length,
            None => unreachable!(),
        };

        let header = Header::builder()
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(SQ0_LN))
            .add_reference_sequence("sq1", Map::<ReferenceSequence>::new(SQ1_LN))
            .build();

        assert_eq!(parse_reference_sequence_id(&header, b"sq0"), Ok(0));
        assert_eq!(parse_reference_sequence_id(&header, b"sq1"), Ok(1));

        assert_eq!(
            parse_reference_sequence_id(&header, b"*"),
            Err(ParseError::MissingReferenceSequenceDictionaryEntry(
                Vec::from("*")
            ))
        );

        assert_eq!(
            parse_reference_sequence_id(&header, b"sq2"),
            Err(ParseError::MissingReferenceSequenceDictionaryEntry(
                Vec::from("sq2")
            ))
        );
    }
}
