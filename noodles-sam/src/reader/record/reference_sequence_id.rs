use std::{io, str};

use crate::Header;

pub(super) fn parse_reference_sequence_id(header: &Header, src: &[u8]) -> io::Result<usize> {
    str::from_utf8(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|s| {
            header.reference_sequences().get_index_of(s).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid reference sequence name",
                )
            })
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

        assert!(matches!(
            parse_reference_sequence_id(&header, b"*"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        assert!(matches!(
            parse_reference_sequence_id(&header, b"sq2"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
