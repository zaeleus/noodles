use noodles_sam::alignment::record::Flags;

use super::num::write_u16_le;

pub(super) fn write_flags(dst: &mut Vec<u8>, flags: Flags) {
    write_u16_le(dst, u16::from(flags));
}
