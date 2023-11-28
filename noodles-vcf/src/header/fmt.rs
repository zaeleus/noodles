use std::fmt;

use super::record::value::{map::Other, Map};

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

pub(super) fn write_meta_record(
    f: &mut fmt::Formatter<'_>,
    id: &str,
    map: &Map<Other>,
) -> fmt::Result {
    use super::record::PREFIX;

    const NUMBER: &str = "Number";
    const TYPE: &str = "Type";
    const VALUES: &str = "Values";

    write!(f, "{}META=<{}={}", PREFIX, map.id_tag(), id)?;

    for (k, v) in map.other_fields() {
        match k.as_ref() {
            NUMBER | TYPE | VALUES => write!(f, ",{k}={v}")?,
            _ => {
                write!(f, ",{k}=")?;
                write_escaped_string(f, v)?;
            }
        }
    }

    write!(f, ">")?;

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

    #[test]
    fn test_write_meta_record() -> Result<(), Box<dyn std::error::Error>> {
        struct Record(&'static str, Map<Other>);

        impl fmt::Display for Record {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write_meta_record(f, self.0, &self.1)
            }
        }

        let record = Record(
            "Assay",
            Map::<Other>::builder()
                .insert("Type".parse()?, "String")
                .insert("Number".parse()?, ".")
                .insert("Values".parse()?, "[WholeGenome, Exome]")
                .build()?,
        );

        let expected = r#"##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>"#;

        assert_eq!(record.to_string(), expected);

        Ok(())
    }
}
