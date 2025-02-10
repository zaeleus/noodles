use std::io::{self, Write};

use crate::MAGIC_NUMBER;

pub(super) fn write_magic_number<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&MAGIC_NUMBER)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_magic_number() -> io::Result<()> {
        let mut buf = Vec::new();
        write_magic_number(&mut buf)?;
        assert_eq!(buf, [b'C', b'R', b'A', b'M']);
        Ok(())
    }
}
