use std::io::{self, Write};

use crate::record::data::field::value::Subtype;

pub fn write_subtype<W>(writer: &mut W, subtype: Subtype) -> io::Result<()>
where
    W: Write,
{
    let n = u8::from(subtype);
    writer.write_all(&[n])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_subtype() -> io::Result<()> {
        let mut buf = Vec::new();
        write_subtype(&mut buf, Subtype::Int32)?;
        assert_eq!(buf, b"i");
        Ok(())
    }
}
