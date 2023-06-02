use bytes::BufMut;
use noodles_sam::{self as sam, record::MappingQuality};

pub fn put_mapping_quality<B>(dst: &mut B, mapping_quality: Option<MappingQuality>)
where
    B: BufMut,
{
    use sam::record::mapping_quality::MISSING;
    let mapq = mapping_quality.map(u8::from).unwrap_or(MISSING);
    dst.put_u8(mapq);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_mapping_quality() {
        fn t(buf: &mut Vec<u8>, mapping_quality: Option<MappingQuality>, expected: &[u8]) {
            buf.clear();
            put_mapping_quality(buf, mapping_quality);
            assert_eq!(buf, expected);
        }

        let mut buf = Vec::new();

        t(&mut buf, None, &[0xff]);
        t(&mut buf, MappingQuality::new(8), &[0x08]);
    }
}
