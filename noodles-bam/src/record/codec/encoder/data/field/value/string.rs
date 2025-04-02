use crate::record::codec::encoder::num::write_u8;

pub(super) fn write_string(dst: &mut Vec<u8>, buf: &[u8]) {
    const NUL: u8 = 0x00;

    dst.extend(buf);
    write_u8(dst, NUL);
}
