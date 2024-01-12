use std::io::{self, Write};

use crate::{alignment::record_buf::Flags, io::writer::num};

pub(super) fn write_flags<W>(writer: &mut W, flags: Flags) -> io::Result<()>
where
    W: Write,
{
    num::write_u16(writer, u16::from(flags))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_flags() -> io::Result<()> {
        let mut buf = Vec::new();
        write_flags(&mut buf, Flags::UNMAPPED)?;
        assert_eq!(buf, b"4");
        Ok(())
    }
}
