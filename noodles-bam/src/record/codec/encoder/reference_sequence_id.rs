use std::{error, fmt};

use noodles_sam as sam;

use super::num::write_i32_le;

const MAX_REFERENCE_SEQUENCE_ID: usize = i32::MAX as usize;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EncodeError {
    OutOfRange(usize),
    MissingEntry { actual: usize, expected: usize },
}

impl error::Error for EncodeError {}

impl fmt::Display for EncodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::OutOfRange(actual) => {
                write!(
                    f,
                    "out of range: expected <= {MAX_REFERENCE_SEQUENCE_ID}, got {actual}"
                )
            }
            Self::MissingEntry { actual, expected } => {
                write!(f, "missing entry: expected < {expected}, got {actual}")
            }
        }
    }
}

pub(super) fn write_reference_sequence_id(
    dst: &mut Vec<u8>,
    header: &sam::Header,
    reference_sequence_id: Option<usize>,
) -> Result<(), EncodeError> {
    const UNMAPPED: i32 = -1;

    let n = if let Some(id) = reference_sequence_id {
        if id < header.reference_sequences().len() {
            i32::try_from(id).map_err(|_| EncodeError::OutOfRange(id))?
        } else {
            return Err(EncodeError::MissingEntry {
                actual: id,
                expected: header.reference_sequences().len(),
            });
        }
    } else {
        UNMAPPED
    };

    write_i32_le(dst, n);

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use super::*;

    #[test]
    fn test_write_reference_sequence_id() -> Result<(), EncodeError> {
        use sam::header::record::value::{Map, map::ReferenceSequence};

        fn t(
            buf: &mut Vec<u8>,
            header: &sam::Header,
            reference_sequence_id: Option<usize>,
            expected: &[u8],
        ) -> Result<(), EncodeError> {
            buf.clear();
            write_reference_sequence_id(buf, header, reference_sequence_id)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let header = sam::Header::default();
        let reference_sequence_id = None;
        t(
            &mut buf,
            &header,
            reference_sequence_id,
            &[0xff, 0xff, 0xff, 0xff],
        )?;

        let header = sam::Header::builder()
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
            )
            .build();
        let reference_sequence_id = Some(0);
        t(
            &mut buf,
            &header,
            reference_sequence_id,
            &[0x00, 0x00, 0x00, 0x00],
        )?;

        buf.clear();
        let header = sam::Header::default();
        let reference_sequence_id = Some(0);
        assert_eq!(
            write_reference_sequence_id(&mut buf, &header, reference_sequence_id),
            Err(EncodeError::MissingEntry {
                actual: 0,
                expected: 0
            })
        );

        Ok(())
    }
}
