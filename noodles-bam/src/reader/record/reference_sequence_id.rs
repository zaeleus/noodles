use std::{io, mem};

use bytes::Buf;

pub(crate) fn get_reference_sequence_id<B>(src: &mut B, n_ref: usize) -> io::Result<Option<usize>>
where
    B: Buf,
{
    const UNMAPPED: i32 = -1;

    if src.remaining() < mem::size_of::<i32>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    match src.get_i32_le() {
        UNMAPPED => Ok(None),
        n => usize::try_from(n)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|m| {
                if m < n_ref {
                    Ok(Some(m))
                } else {
                    Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "invalid reference sequence ID: expected < {}, got {}",
                            n_ref, m
                        ),
                    ))
                }
            }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_reference_sequence_id() -> io::Result<()> {
        let data = (-1i32).to_le_bytes();
        let mut src = &data[..];
        assert!(get_reference_sequence_id(&mut src, 1)?.is_none());

        let data = 0i32.to_le_bytes();
        let mut src = &data[..];
        assert_eq!(get_reference_sequence_id(&mut src, 1)?, Some(0));

        let data = [];
        let mut src = &data[..];
        assert!(matches!(
            get_reference_sequence_id(&mut src, 1),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = (-2i32).to_le_bytes();
        let mut src = &data[..];
        assert!(matches!(
            get_reference_sequence_id(&mut src, 1),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        let data = 0i32.to_le_bytes();
        let mut src = &data[..];
        assert!(matches!(
            get_reference_sequence_id(&mut src, 0),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
