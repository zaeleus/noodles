use bytes::BufMut;
use noodles_sam::alignment::record::Flags;

pub(super) fn put_flags<B>(dst: &mut B, flags: Flags)
where
    B: BufMut,
{
    dst.put_u16_le(u16::from(flags));
}
