use std::io;

use super::read_value;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Type {
    Int8(usize),
    Int16(usize),
    Int32(usize),
    Float(usize),
    String(usize),
}

const MAX_TYPE_LEN: usize = 0x0f;

pub(crate) fn read_type(src: &mut &[u8]) -> io::Result<Option<Type>> {
    let encoding = get_u8(src)?;
    let mut len = usize::from(encoding >> 4);

    if len == MAX_TYPE_LEN {
        let value = read_value(src)?;

        len = match value.and_then(|v| v.as_int()) {
            Some(n) => {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
            }
            None => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid length value",
                ));
            }
        };
    }

    match encoding & 0x0f {
        0 => Ok(None),
        1 => Ok(Some(Type::Int8(len))),
        2 => Ok(Some(Type::Int16(len))),
        3 => Ok(Some(Type::Int32(len))),
        5 => Ok(Some(Type::Float(len))),
        7 => Ok(Some(Type::String(len))),
        _ => Err(io::Error::new(io::ErrorKind::InvalidData, "invalid type")),
    }
}

fn get_u8(src: &mut &[u8]) -> io::Result<u8> {
    if let Some((b, rest)) = src.split_first() {
        *src = rest;
        Ok(*b)
    } else {
        Err(io::Error::from(io::ErrorKind::UnexpectedEof))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_type() -> io::Result<()> {
        fn t(mut src: &[u8], expected: Option<Type>) -> io::Result<()> {
            assert_eq!(read_type(&mut src)?, expected);
            Ok(())
        }

        t(&[0x00], None)?;
        t(&[0x10], None)?;
        t(&[0x11], Some(Type::Int8(1)))?;
        t(&[0x12], Some(Type::Int16(1)))?;
        t(&[0x13], Some(Type::Int32(1)))?;
        t(&[0x15], Some(Type::Float(1)))?;
        t(&[0x17], Some(Type::String(1)))?;

        t(
            &[
                0xf1, // (len >= 15, Type::Int8)
                0x11, // Some(Type::Int8(1))
                0x15, // Some(Value::Int8(21))
            ],
            Some(Type::Int8(21)),
        )?;

        let mut src = &[][..];
        assert!(matches!(read_type(&mut src), Err(e) if e.kind() == io::ErrorKind::UnexpectedEof));

        let mut src = &[0x14][..];
        assert!(matches!(read_type(&mut src), Err(e) if e.kind() == io::ErrorKind::InvalidData));

        Ok(())
    }
}
