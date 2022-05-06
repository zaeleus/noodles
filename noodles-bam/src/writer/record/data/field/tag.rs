use bytes::BufMut;
use noodles_sam::record::data::field::Tag;

pub fn put_tag<B>(dst: &mut B, tag: Tag)
where
    B: BufMut,
{
    dst.put(&tag.as_ref()[..]);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_tag() {
        let mut buf = Vec::new();
        put_tag(&mut buf, Tag::AlignmentHitCount);
        assert_eq!(buf, [b'N', b'H']);
    }
}
