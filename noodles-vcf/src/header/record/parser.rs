use nom::{
    branch::alt,
    bytes::complete::{escaped_transform, tag, take_till, take_until},
    character::complete::{alpha1, none_of},
    combinator::{map, rest},
    multi::separated_nonempty_list,
    sequence::{delimited, separated_pair},
    IResult,
};

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

fn structure(input: &str) -> IResult<&str, Vec<(String, String)>> {
    delimited(tag("<"), separated_nonempty_list(tag(","), field), tag(">"))(input)
}

pub fn parse(input: &str) -> IResult<&str, (String, Vec<(String, String)>)> {
    let (input, _) = tag("##")(input)?;
    let (input, key) = take_until("=")(input)?;
    let (input, _) = tag("=")(input)?;

    let (input, values) = if input.starts_with("<") {
        structure(input)?
    } else {
        map(rest, |s: &str| vec![(key.into(), s.into())])(input)?
    };

    Ok((input, (key.into(), values)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() -> Result<(), Box<dyn std::error::Error>> {
        let (_, (key, values)) = parse("##fileformat=VCFv4.3")?;
        assert_eq!(key, "fileformat");
        assert_eq!(
            values,
            vec![(String::from("fileformat"), String::from("VCFv4.3"))]
        );

        let (_, (key, values)) = parse("##fileDate=20200502")?;
        assert_eq!(key, "fileDate");
        assert_eq!(
            values,
            vec![(String::from("fileDate"), String::from("20200502"))]
        );

        let (_, (key, values)) = parse("##reference=file:///tmp/ref.fasta")?;
        assert_eq!(key, "reference");
        assert_eq!(
            values,
            vec![(
                String::from("reference"),
                String::from("file:///tmp/ref.fasta")
            )]
        );

        let (_, (key, values)) = parse(
            r#"##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">"#,
        )?;
        assert_eq!(key, "INFO");
        assert_eq!(
            values,
            vec![
                (String::from("ID"), String::from("NS")),
                (String::from("Number"), String::from("1")),
                (String::from("Type"), String::from("Integer")),
                (
                    String::from("Description"),
                    String::from("Number of samples with data")
                ),
            ]
        );

        assert!(parse("").is_err());
        assert!(parse("fileformat=VCFv4.3").is_err());
        assert!(parse("#fileformat=VCFv4.3").is_err());

        Ok(())
    }
}
