use bytes::BufMut;
use noodles_sam as sam;

pub fn put_read_name<B>(dst: &mut B, read_name: Option<&sam::record::ReadName>)
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
