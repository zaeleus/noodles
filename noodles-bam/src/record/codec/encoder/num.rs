pub(super) fn write_i8(dst: &mut Vec<u8>, n: i8) {
    write_u8(dst, n as u8);
}

pub(super) fn write_u8(dst: &mut Vec<u8>, n: u8) {
    dst.push(n);
}

pub(super) fn write_i16_le(dst: &mut Vec<u8>, n: i16) {
    dst.extend(n.to_le_bytes());
}

pub(super) fn write_u16_le(dst: &mut Vec<u8>, n: u16) {
    dst.extend(n.to_le_bytes());
}

pub(super) fn write_i32_le(dst: &mut Vec<u8>, n: i32) {
    dst.extend(n.to_le_bytes());
}

pub(super) fn write_u32_le(dst: &mut Vec<u8>, n: u32) {
    dst.extend(n.to_le_bytes());
}

pub(super) fn write_f32_le(dst: &mut Vec<u8>, n: f32) {
    dst.extend(n.to_le_bytes());
}
