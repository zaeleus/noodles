use nom::{
    branch::alt,
    bytes::complete::{escaped_transform, tag, take_till, take_until},
    character::complete::{alpha1, alphanumeric1, none_of},
    combinator::{map, opt, rest},
    multi::separated_nonempty_list,
    sequence::{delimited, separated_pair},
    IResult,
};

use super::Value;

pub enum Line {
    Meta(String, Value),
    Header(Vec<String>),
}

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

fn meta(input: &str) -> IResult<&str, Line> {
    let (input, key) = delimited(tag("##"), take_until("="), tag("="))(input)?;
    let (input, value) = alt((structure, map(rest, |s: &str| Value::String(s.into()))))(input)?;
    Ok((input, Line::Meta(key.into(), value)))
}

fn header(input: &str) -> IResult<&str, Line> {
    let (input, _) = tag("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")(input)?;
    let (input, format) = opt(tag("\tFORMAT\t"))(input)?;

    if format.is_some() {
        let (input, values) = separated_nonempty_list(tag("\t"), alphanumeric1)(input)?;
        let sample_names = values.into_iter().map(|s| s.into()).collect();
        Ok((input, Line::Header(sample_names)))
    } else {
        Ok((input, Line::Header(Vec::new())))
    }
}

pub fn parse(input: &str) -> IResult<&str, Line> {
    alt((meta, header))(input)
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
    fn test_parse_with_meta() -> Result<(), Box<dyn std::error::Error>> {
        let (_, line) = parse("##fileformat=VCFv4.3")?;

        match line {
            Line::Meta(key, value) => {
                assert_eq!(key, "fileformat");
                assert_eq!(value, Value::String(String::from("VCFv4.3")));
            }
            Line::Header(_) => panic!(),
        }

        let (_, line) = parse("##fileDate=20200502")?;

        match line {
            Line::Meta(key, value) => {
                assert_eq!(key, "fileDate");
                assert_eq!(value, Value::String(String::from("20200502")));
            }
            Line::Header(_) => panic!(),
        }

        let (_, line) = parse("##reference=file:///tmp/ref.fasta")?;

        match line {
            Line::Meta(key, value) => {
                assert_eq!(key, "reference");
                assert_eq!(value, Value::String(String::from("file:///tmp/ref.fasta")));
            }
            Line::Header(_) => panic!(),
        }

        let (_, line) = parse(
            r#"##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">"#,
        )?;

        match line {
            Line::Meta(key, value) => {
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
            }
            Line::Header(_) => panic!(),
        }

        Ok(())
    }

    #[test]
    fn test_parse_with_header() -> Result<(), Box<dyn std::error::Error>> {
        let (_, line) = parse("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;

        match line {
            Line::Header(sample_names) => {
                assert!(sample_names.is_empty());
            }
            Line::Meta(..) => panic!(),
        }

        let (_, line) = parse("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample0")?;

        match line {
            Line::Header(sample_names) => {
                assert_eq!(sample_names.len(), 1);
                assert_eq!(&sample_names[0], "sample0");
            }
            Line::Meta(..) => panic!(),
        }

        let (_, line) =
            parse("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample0\tsample1")?;

        match line {
            Line::Header(sample_names) => {
                assert_eq!(sample_names.len(), 2);
                assert_eq!(&sample_names[0], "sample0");
                assert_eq!(&sample_names[1], "sample1");
            }
            Line::Meta(..) => panic!(),
        }

        Ok(())
    }
}
