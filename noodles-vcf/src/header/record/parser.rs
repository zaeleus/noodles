use nom::{
    branch::alt,
    bytes::complete::{escaped_transform, tag, take_till, take_until},
    character::complete::{alpha1, none_of},
    combinator::{map, rest},
    multi::separated_nonempty_list,
    sequence::{delimited, separated_pair},
    IResult,
};

use super::{Value, PREFIX};

fn string(input: &str) -> IResult<&str, String> {
    delimited(
        tag("\""),
        escaped_transform(none_of("\\\""), '\\', alt((tag("\\"), tag("\"")))),
        tag("\""),
    )(input)
}

fn value(input: &str) -> IResult<&str, String> {
    map(take_till(|c| c == ',' || c == '>'), |s: &str| s.into())(input)
}

fn field(input: &str) -> IResult<&str, (String, String)> {
    map(
        separated_pair(alpha1, tag("="), alt((string, value))),
        |(k, v): (&str, String)| (k.into(), v),
    )(input)
}

fn structure(input: &str) -> IResult<&str, Value> {
    map(
        delimited(tag("<"), separated_nonempty_list(tag(","), field), tag(">")),
        Value::Struct,
    )(input)
}

fn meta(input: &str) -> IResult<&str, (String, Value)> {
    let (input, key) = delimited(tag(PREFIX), take_until("="), tag("="))(input)?;
    let (input, value) = alt((structure, map(rest, |s: &str| Value::String(s.into()))))(input)?;
    Ok((input, (key.into(), value)))
}

pub fn parse(input: &str) -> IResult<&str, (String, Value)> {
    meta(input)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() {
        assert!(parse("").is_err());
        assert!(parse("fileformat=VCFv4.3").is_err());
        assert!(parse("#fileformat=VCFv4.3").is_err());
    }

    #[test]
    fn test_parse_with_meta_string_value() -> Result<(), Box<dyn std::error::Error>> {
        let (_, (key, value)) = parse("##fileformat=VCFv4.3")?;
        assert_eq!(key, "fileformat");
        assert_eq!(value, Value::String(String::from("VCFv4.3")));

        let (_, (key, value)) = parse("##fileDate=20200502")?;
        assert_eq!(key, "fileDate");
        assert_eq!(value, Value::String(String::from("20200502")));

        let (_, (key, value)) = parse("##reference=file:///tmp/ref.fasta")?;
        assert_eq!(key, "reference");
        assert_eq!(value, Value::String(String::from("file:///tmp/ref.fasta")));

        Ok(())
    }

    #[test]
    fn test_parse_with_meta_struct_value() -> Result<(), Box<dyn std::error::Error>> {
        let (_, (key, value)) = parse(
            r#"##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">"#,
        )?;

        assert_eq!(key, "ALT");
        assert_eq!(
            value,
            Value::Struct(vec![
                (String::from("ID"), String::from("NON_REF")),
                (
                    String::from("Description"),
                    String::from("Represents any possible alternative allele at this location")
                ),
            ])
        );

        let (_, (key, value)) = parse(
            r#"##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">"#,
        )?;

        assert_eq!(key, "INFO");
        assert_eq!(
            value,
            Value::Struct(vec![
                (String::from("ID"), String::from("NS")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("Integer")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data")
                ),
            ])
        );

        Ok(())
    }
}
