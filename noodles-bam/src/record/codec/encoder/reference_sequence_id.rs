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
