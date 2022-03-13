use std::io::{self, Write};

use byteorder::WriteBytesExt;
use noodles_sam::record::data::field::value::Subtype;

pub fn write_subtype<W>(writer: &mut W, subtype: Subtype) -> io::Result<()>
where
    W: Write,
{
    writer.write_u8(u8::from(subtype))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_type() -> io::Result<()> {
        let mut buf = Vec::new();
        write_subtype(&mut buf, Subtype::Int32)?;
        assert_eq!(buf, [b'i']);
        Ok(())
    }
}
