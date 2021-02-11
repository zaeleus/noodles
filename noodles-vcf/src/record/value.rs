use std::{borrow::Cow, num, str};

use percent_encoding::percent_decode_str;

/// Parses a single-precision floating-point.
///
/// This is extended to support case-insensitive values of (+/-)infinity and NaN.
///
/// See § 1.3 Data types (2020-06-25): "Float (32-bit IEEE-754, formatted to match [...]
/// ^[-+]?(INF|INFINITY|NAN)$ case insensitively)".
pub(crate) fn parse_f32_case_insensitive_extended(s: &str) -> Result<f32, num::ParseFloatError> {
    match s.to_lowercase().as_ref() {
        "-infinity" => Ok(f32::NEG_INFINITY),
        "infinity" | "+infinity" => Ok(f32::INFINITY),
        "nan" | "-nan" | "+nan" => Ok(f32::NAN),
        t => t.parse(),
    }
}

pub(crate) fn percent_decode(s: &str) -> Result<Cow<str>, str::Utf8Error> {
    percent_decode_str(s).decode_utf8()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_f32_case_insensitive_extended() -> Result<(), num::ParseFloatError> {
        assert_eq!(
            parse_f32_case_insensitive_extended("-inf"),
            Ok(f32::NEG_INFINITY)
        );
        assert_eq!(
            parse_f32_case_insensitive_extended("-infinity"),
            Ok(f32::NEG_INFINITY)
        );
        assert_eq!(
            parse_f32_case_insensitive_extended("+inf"),
            Ok(f32::INFINITY)
        );
        assert_eq!(
            parse_f32_case_insensitive_extended("+infinity"),
            Ok(f32::INFINITY)
        );
        assert_eq!(
            parse_f32_case_insensitive_extended("infinity"),
            Ok(f32::INFINITY)
        );
        assert_eq!(
            parse_f32_case_insensitive_extended("INFINITY"),
            Ok(f32::INFINITY)
        );

        assert!(parse_f32_case_insensitive_extended("-NaN")?.is_nan());
        assert!(parse_f32_case_insensitive_extended("+NaN")?.is_nan());
        assert!(parse_f32_case_insensitive_extended("NaN")?.is_nan());
        assert!(parse_f32_case_insensitive_extended("nan")?.is_nan());

        assert_eq!(parse_f32_case_insensitive_extended("0.0"), Ok(0.0));

        Ok(())
    }
}
