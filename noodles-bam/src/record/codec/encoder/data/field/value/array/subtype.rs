use bytes::BufMut;
use noodles_sam::alignment::record::data::field::value::array::Subtype;

pub fn put_subtype<B>(dst: &mut B, subtype: Subtype)
where
    B: BufMut,
{
    let n = match subtype {
        Subtype::Int8 => b'c',
        Subtype::UInt8 => b'C',
        Subtype::Int16 => b's',
        Subtype::UInt16 => b'S',
        Subtype::Int32 => b'i',
        Subtype::UInt32 => b'I',
        Subtype::Float => b'f',
    };

    dst.put_u8(n);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_subtype() {
        fn t(buf: &mut Vec<u8>, subtype: Subtype, expected: &[u8]) {
            buf.clear();
            put_subtype(buf, subtype);
            assert_eq!(buf, expected);
        }

        let mut buf = Vec::new();

        t(&mut buf, Subtype::Int8, &[b'c']);
        t(&mut buf, Subtype::UInt8, &[b'C']);
        t(&mut buf, Subtype::Int16, &[b's']);
        t(&mut buf, Subtype::UInt16, &[b'S']);
        t(&mut buf, Subtype::Int32, &[b'i']);
        t(&mut buf, Subtype::UInt32, &[b'I']);
        t(&mut buf, Subtype::Float, &[b'f']);
    }
}
