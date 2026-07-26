use std::borrow::Cow;

use bstr::{BStr, ByteSlice};

pub(super) fn decode(s: &[u8]) -> Cow<'_, BStr> {
    match Cow::from(percent_encoding::percent_decode(s)) {
        Cow::Borrowed(buf) => Cow::Borrowed(buf.as_bstr()),
        Cow::Owned(buf) => Cow::Owned(buf.into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode() {
        assert_eq!(decode(b"sq0"), BStr::new("sq0"));
        assert_eq!(decode(b"sq%200"), BStr::new("sq 0"));
        assert_eq!(decode(b"100%25"), BStr::new("100%"));
        assert_eq!(decode(b"a%09b"), BStr::new("a\tb"));

        // Invalid escapes are passed through.
        assert_eq!(decode(b"50%"), BStr::new("50%"));
        assert_eq!(decode(b"%zz"), BStr::new("%zz"));
    }
}
