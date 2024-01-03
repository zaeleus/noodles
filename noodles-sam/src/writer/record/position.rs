use std::io::{self, Write};

use crate::{alignment::record::Position, writer::num};

// § 1.4 "The alignment section: mandatory fields" (2021-06-03): `[0, 2^31-1]`.
const MAX_POSITION_VALUE: usize = (1 << 31) - 1;

pub fn write_position<W>(writer: &mut W, position: Option<&dyn Position>) -> io::Result<()>
where
    W: Write,
{
    let pos = position
        .map(|position| position.try_to_usize())
        .transpose()?
        .unwrap_or_default();

    if pos <= MAX_POSITION_VALUE {
        num::write_usize(writer, pos)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid position: expected position to be <= {MAX_POSITION_VALUE}, got {pos}"),
        ))
    }
}

#[cfg(test)]
mod tests {
    use noodles_core as core;

    use super::*;

    #[test]
    fn test_write_position() -> io::Result<()> {
        fn t(
            buf: &mut Vec<u8>,
            position_buf: Option<core::Position>,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();

            let position = position_buf
                .as_ref()
                .map(|position| position as &dyn Position);

            write_position(buf, position)?;

            assert_eq!(buf, expected);

            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, b"0")?;
        t(&mut buf, core::Position::new(13), b"13")?;

        buf.clear();
        let position_buf = core::Position::new(MAX_POSITION_VALUE + 1);
        let position = position_buf
            .as_ref()
            .map(|position| position as &dyn Position);
        assert!(matches!(
            write_position(&mut buf, position),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }
}
