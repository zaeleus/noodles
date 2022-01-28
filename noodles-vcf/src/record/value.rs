use std::{borrow::Cow, num, str};

use percent_encoding::percent_decode_str;

/// Parses a single-precision floating-point.
pub(crate) fn parse_f32(s: &str) -> Result<f32, num::ParseFloatError> {
    s.parse()
}

pub(crate) fn percent_decode(s: &str) -> Result<Cow<'_, str>, str::Utf8Error> {
    percent_decode_str(s).decode_utf8()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_f32() -> Result<(), num::ParseFloatError> {
        assert_eq!(parse_f32("-inf"), Ok(f32::NEG_INFINITY));
        assert_eq!(parse_f32("-infinity"), Ok(f32::NEG_INFINITY));
        assert_eq!(parse_f32("+inf"), Ok(f32::INFINITY));
        assert_eq!(parse_f32("+infinity"), Ok(f32::INFINITY));
        assert_eq!(parse_f32("infinity"), Ok(f32::INFINITY));
        assert_eq!(parse_f32("INFINITY"), Ok(f32::INFINITY));

        assert!(parse_f32("-NaN")?.is_nan());
        assert!(parse_f32("+NaN")?.is_nan());
        assert!(parse_f32("NaN")?.is_nan());
        assert!(parse_f32("nan")?.is_nan());

        assert_eq!(parse_f32("0.0"), Ok(0.0));

        Ok(())
    }

    #[test]
    fn test_percent_decode() -> Result<(), str::Utf8Error> {
        assert_eq!(percent_decode("noodles")?, "noodles");
        assert_eq!(percent_decode("noodles%3Dvcf")?, "noodles=vcf");
        Ok(())
    }
}
