use std::{borrow::Cow, io};

use super::percent_decode;

pub(super) fn parse_tag(s: &str) -> io::Result<Cow<'_, str>> {
    percent_decode(s).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tag() -> io::Result<()> {
        assert_eq!(parse_tag("ID")?, Cow::from("ID"));
        assert_eq!(parse_tag("%25s")?, Cow::from("%s"));
        Ok(())
    }
}
