use std::fmt;

/// Formats a string as an escaped string.
///
/// This writes double quote delimiters and escapes two characters: the double quote (`"`) and
/// backslash (`\`).
pub(crate) fn write_escaped_string(f: &mut fmt::Formatter<'_>, s: &str) -> fmt::Result {
    f.write_str("\"")?;

    for c in s.chars() {
        if matches!(c, '"' | '\\') {
            f.write_str("\\")?;
        }

        write!(f, "{c}")?;
    }

    f.write_str("\"")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    struct EscapedStringFormat(&'static str);

    impl fmt::Display for EscapedStringFormat {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write_escaped_string(f, self.0)
        }
    }

    #[test]
    fn test_write_escaped_string() {
        assert_eq!(
            EscapedStringFormat(r#"noodles"#).to_string(),
            r#""noodles""#
        );

        assert_eq!(
            EscapedStringFormat("noodles=🍜").to_string(),
            r#""noodles=🍜""#
        );

        assert_eq!(
            EscapedStringFormat(r#"noodles-"vcf""#).to_string(),
            r#""noodles-\"vcf\"""#
        );

        assert_eq!(
            EscapedStringFormat(r"noodles\vcf").to_string(),
            r#""noodles\\vcf""#
        );
    }
}
