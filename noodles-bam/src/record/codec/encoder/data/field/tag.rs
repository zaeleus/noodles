use bytes::BufMut;
use noodles_sam::alignment::record::data::field::Tag;

pub fn put_tag<B>(dst: &mut B, tag: Tag)
where
    B: BufMut,
{
    dst.put_slice(tag.as_ref());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_tag() {
        let mut buf = Vec::new();
        put_tag(&mut buf, Tag::ALIGNMENT_HIT_COUNT);
        assert_eq!(buf, [b'N', b'H']);
    }
}
