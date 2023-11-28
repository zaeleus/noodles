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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_position() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, position: Option<Position>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            put_position(buf, position)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, &[0xff, 0xff, 0xff, 0xff])?;
        t(&mut buf, Some(Position::MIN), &[0x00, 0x00, 0x00, 0x00])?;
        t(
            &mut buf,
            Position::try_from(8).map(Some)?,
            &[0x07, 0x00, 0x00, 0x00],
        )?;

        Ok(())
    }
}
