use std::io;

use bytes::BufMut;
use noodles_core::Position;

pub(super) fn put_position<B>(dst: &mut B, position: Option<Position>) -> io::Result<()>
where
    B: BufMut,
{
    const MISSING: i32 = -1;

    let pos = if let Some(position) = position {
        i32::try_from(usize::from(position) - 1)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
    } else {
        MISSING
    };

    dst.put_i32_le(pos);

    Ok(())
}
