use std::borrow::Cow;

use bstr::BStr;

use super::percent_decode;

pub(super) fn parse_tag(src: &[u8]) -> Cow<'_, BStr> {
    percent_decode(src)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tag() {
        assert_eq!(parse_tag(b"ID"), Cow::from(BStr::new("ID")));
        assert_eq!(parse_tag(b"%25s"), Cow::from(BStr::new("%s")));
    }
}
