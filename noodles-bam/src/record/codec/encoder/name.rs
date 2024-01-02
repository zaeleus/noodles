use bytes::BufMut;
use noodles_sam::{self as sam, record::Name};

pub fn put_name<B>(dst: &mut B, name: Option<&Name>)
where
    B: BufMut,
{
    use sam::record::name::MISSING;

    const NUL: u8 = 0x00;

    if let Some(name) = name {
        dst.put(name.as_ref());
    } else {
        dst.put(MISSING);
    }

    dst.put_u8(NUL);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_name() -> Result<(), sam::record::name::ParseError> {
        fn t(buf: &mut Vec<u8>, name: Option<&Name>, expected: &[u8]) {
            buf.clear();
            put_name(buf, name);
            assert_eq!(buf, expected);
        }

        let mut buf = Vec::new();

        t(&mut buf, None, &[b'*', 0x00]);
        t(&mut buf, Some(&"r".parse()?), &[b'r', 0x00]);
        t(&mut buf, Some(&"r0".parse()?), &[b'r', b'0', 0x00]);

        Ok(())
    }
}
