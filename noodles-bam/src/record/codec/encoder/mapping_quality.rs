use noodles_sam::alignment::record::MappingQuality;

use super::num::write_u8;

pub(super) fn write_mapping_quality(dst: &mut Vec<u8>, mapping_quality: Option<MappingQuality>) {
    const MISSING: u8 = 0xff;
    let n = mapping_quality.map(u8::from).unwrap_or(MISSING);
    write_u8(dst, n);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_mapping_quality() {
        fn t(buf: &mut Vec<u8>, mapping_quality: Option<MappingQuality>, expected: &[u8]) {
            buf.clear();
            write_mapping_quality(buf, mapping_quality);
            assert_eq!(buf, expected);
        }

        let mut buf = Vec::new();

        t(&mut buf, None, &[0xff]);
        t(&mut buf, MappingQuality::new(8), &[0x08]);
    }
}
