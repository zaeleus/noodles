use std::io::{self, Write};

use noodles_core::Position;

use crate::io::writer::num::write_i32_le;

pub(super) fn write_position<W>(writer: &mut W, position: Option<Position>) -> io::Result<()>
where
    W: Write,
{
    const TELOMERE_START: i32 = -1;

    let pos = match position {
        Some(position) => i32::try_from(usize::from(position))
            .map(|n| n - 1)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?,
        None => TELOMERE_START,
    };

    write_i32_le(writer, pos)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_pos() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, position: Option<Position>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_position(buf, position)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, &[0xff, 0xff, 0xff, 0xff])?;
        t(&mut buf, Some(Position::MIN), &[0x00, 0x00, 0x00, 0x00])?;
        t(
            &mut buf,
            Some(Position::try_from((1 << 31) - 1)?),
            &[0xfe, 0xff, 0xff, 0x7f],
        )?;

        buf.clear();
        assert!(matches!(
            write_position(&mut buf, Some(Position::try_from(1 << 32)?)),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
