use std::{io, num::NonZeroUsize};

use bytes::Buf;
use noodles_sam as sam;

pub(super) fn get_read_name<B>(
    buf: &mut B,
    read_name: &mut Option<sam::record::ReadName>,
    l_read_name: NonZeroUsize,
) -> io::Result<()>
where
    B: Buf,
{
    const MISSING: u8 = b'*';

    let len = usize::from(l_read_name);

    if buf.remaining() < len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    *read_name = if len == 2 && buf.chunk()[0] == MISSING {
        buf.advance(2);
        None
    } else {
        let mut read_name_buf = read_name.take().map(Vec::from).unwrap_or_default();

        // SAFETY: len is guaranteed to be > 0.
        read_name_buf.resize(len - 1, Default::default());
        buf.copy_to_slice(&mut read_name_buf);

        // Discard the NUL terminator.
        buf.advance(1);

        sam::record::ReadName::try_from(read_name_buf)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
    };

    Ok(())
}
