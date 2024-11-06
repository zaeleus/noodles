use std::io;

use lexical_core::FromLexical;

use crate::alignment::record::data::field::Value;

pub(super) fn parse_integer_value<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    parse_int(src)
        .map(Value::Int32)
        .or_else(|e| match e {
            lexical_core::Error::Overflow(_) => parse_int(src).map(Value::UInt32),
            _ => Err(e),
        })
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn parse_int<N>(src: &mut &[u8]) -> lexical_core::Result<N>
where
    N: FromLexical,
{
    let (n, i) = lexical_core::parse_partial(src)?;
    *src = &src[i..];
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_integer_value() -> io::Result<()> {
        let mut src = &b"-2147483649"[..]; // -(1 << 31) - 1
        assert!(matches!(
            parse_integer_value(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        let mut src = &b"-2147483648"[..];
        let actual = parse_integer_value(&mut src)?;
        assert!(matches!(actual, Value::Int32(i32::MIN)));

        let mut src = &b"0"[..];
        let actual = parse_integer_value(&mut src)?;
        assert!(matches!(actual, Value::Int32(0)));

        let mut src = &b"2147483647"[..];
        let actual = parse_integer_value(&mut src)?;
        assert!(matches!(actual, Value::Int32(i32::MAX)));

        let mut src = &b"2147483648"[..];
        let actual = parse_integer_value(&mut src)?;
        assert!(matches!(actual, Value::UInt32(n) if n == 1 << 31));

        let mut src = &b"4294967295"[..];
        let actual = parse_integer_value(&mut src)?;
        assert!(matches!(actual, Value::UInt32(u32::MAX)));

        let mut src = &b"4294967296"[..]; // 1 << 32
        assert!(matches!(
            parse_integer_value(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
