use std::io;

use bytes::BufMut;
use noodles_core as core;
use noodles_sam::alignment::record::Position;

pub(super) fn put_position<B, P>(dst: &mut B, position: Option<P>) -> io::Result<()>
where
    B: BufMut,
    P: Position,
{
    const MISSING: i32 = -1;

    let pos = if let Some(position) = position {
        let position = core::Position::try_from(&position as &dyn Position)?;

        i32::try_from(usize::from(position) - 1)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
    } else {
        MISSING
    };

    dst.put_i32_le(pos);

    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_put_position() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            buf: &mut Vec<u8>,
            position: Option<core::Position>,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            put_position(buf, position)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, &[0xff, 0xff, 0xff, 0xff])?;
        t(
            &mut buf,
            Some(core::Position::MIN),
            &[0x00, 0x00, 0x00, 0x00],
        )?;
        t(
            &mut buf,
            core::Position::try_from(8).map(Some)?,
            &[0x07, 0x00, 0x00, 0x00],
        )?;

        Ok(())
    }
}
