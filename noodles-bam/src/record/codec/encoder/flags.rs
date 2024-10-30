use noodles_sam::alignment::record::Flags;

use super::num::write_u16_le;

pub(super) fn write_flags(dst: &mut Vec<u8>, flags: Flags) {
    write_u16_le(dst, u16::from(flags));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_flags() {
        let mut buf = Vec::new();

        buf.clear();
        write_flags(&mut buf, Flags::default());
        assert_eq!(buf, [0x00, 0x00]);

        buf.clear();
        write_flags(&mut buf, Flags::UNMAPPED);
        assert_eq!(buf, [0x04, 0x00]);
    }
}
