use bytes::BufMut;

pub fn put_tag<B>(dst: &mut B, tag: [u8; 2])
where
    B: BufMut,
{
    dst.put(&tag[..]);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_tag() {
        let mut buf = Vec::new();
        put_tag(&mut buf, [b'N', b'H']);
        assert_eq!(buf, [b'N', b'H']);
    }
}
