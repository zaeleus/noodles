use bytes::BufMut;
use noodles_sam::{self as sam, record::ReadName};

pub fn put_read_name<B>(dst: &mut B, read_name: Option<&ReadName>)
where
    B: BufMut,
{
    use sam::record::read_name::MISSING;

    const NUL: u8 = 0x00;

    if let Some(read_name) = read_name {
        dst.put(read_name.as_ref());
    } else {
        dst.put(MISSING);
    }

    dst.put_u8(NUL);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_read_name() -> Result<(), sam::record::read_name::ParseError> {
        fn t(buf: &mut Vec<u8>, read_name: Option<&ReadName>, expected: &[u8]) {
            buf.clear();
            put_read_name(buf, read_name);
            assert_eq!(buf, expected);
        }

        let mut buf = Vec::new();

        t(&mut buf, None, &[b'*', 0x00]);
        t(&mut buf, Some(&"r".parse()?), &[b'r', 0x00]);
        t(&mut buf, Some(&"r0".parse()?), &[b'r', b'0', 0x00]);

        Ok(())
    }
}
