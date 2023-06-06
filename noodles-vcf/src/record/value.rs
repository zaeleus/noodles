use std::{borrow::Cow, str};

use percent_encoding::percent_decode_str;

pub(crate) fn percent_decode(s: &str) -> Result<Cow<'_, str>, str::Utf8Error> {
    percent_decode_str(s).decode_utf8()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_percent_decode() -> Result<(), str::Utf8Error> {
        assert_eq!(percent_decode("noodles")?, "noodles");
        assert_eq!(percent_decode("noodles%3Dvcf")?, "noodles=vcf");
        Ok(())
    }
}
