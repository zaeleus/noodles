use std::io;

use bytes::BufMut;
use noodles_sam as sam;

pub(super) fn put_reference_sequence_id<B>(
    dst: &mut B,
    header: &sam::Header,
    reference_sequence_id: Option<usize>,
) -> io::Result<()>
where
    B: BufMut,
{
    const UNMAPPED: i32 = -1;

    let ref_id = if let Some(id) = reference_sequence_id {
        if id < header.reference_sequences().len() {
            i32::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "invalid reference sequence ID: expected < {}, got {}",
                    header.reference_sequences().len(),
                    id
                ),
            ));
        }
    } else {
        UNMAPPED
    };

    dst.put_i32_le(ref_id);

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use super::*;

    #[test]
    fn test_put_reference_sequence_id() -> io::Result<()> {
        use sam::header::record::value::{map::ReferenceSequence, Map};

        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        fn t(
            buf: &mut Vec<u8>,
            header: &sam::Header,
            reference_sequence_id: Option<usize>,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            put_reference_sequence_id(buf, header, reference_sequence_id)?;
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
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(SQ0_LN))
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
        assert!(matches!(
            put_reference_sequence_id(&mut buf, &header, reference_sequence_id),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
