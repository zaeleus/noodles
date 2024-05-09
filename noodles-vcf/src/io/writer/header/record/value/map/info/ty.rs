use std::io::{self, Write};

use crate::header::record::value::map::info::Type;

pub(super) fn write_type<W>(writer: &mut W, ty: Type) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(ty.as_ref().as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_type() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, ty: Type, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_type(buf, ty)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, Type::Integer, b"Integer")?;
        t(&mut buf, Type::Float, b"Float")?;
        t(&mut buf, Type::Flag, b"Flag")?;
        t(&mut buf, Type::Character, b"Character")?;
        t(&mut buf, Type::String, b"String")?;

        Ok(())
    }
}
