use noodles_sam::alignment::record::data::field::Tag;

pub fn write_tag(dst: &mut Vec<u8>, tag: Tag) {
    dst.extend(tag.as_ref());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_tag() {
        let mut buf = Vec::new();
        write_tag(&mut buf, Tag::ALIGNMENT_HIT_COUNT);
        assert_eq!(buf, [b'N', b'H']);
    }
}
