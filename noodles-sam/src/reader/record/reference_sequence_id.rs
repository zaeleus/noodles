use std::{io, str};

use crate::Header;

pub(super) fn parse_reference_sequence_id(
    header: &Header,
    src: &[u8],
) -> io::Result<Option<usize>> {
    const MISSING: &[u8] = b"*";

    match src {
        MISSING => Ok(None),
        _ => str::from_utf8(src)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                header
                    .reference_sequences()
                    .get_index_of(s)
                    .map(Some)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "invalid reference sequence name",
                        )
                    })
            }),
    }
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

        assert!(parse_reference_sequence_id(&header, b"*")?.is_none());

        assert_eq!(parse_reference_sequence_id(&header, b"sq0")?, Some(0));
        assert_eq!(parse_reference_sequence_id(&header, b"sq1")?, Some(1));

        assert!(matches!(
            parse_reference_sequence_id(&header, b"sq2"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
