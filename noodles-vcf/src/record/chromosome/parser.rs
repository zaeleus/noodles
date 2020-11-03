use nom::{
    branch::alt,
    bytes::complete::{tag, take_until, take_while},
    combinator::{all_consuming, map},
    sequence::delimited,
    IResult,
};

#[derive(Debug, Eq, PartialEq)]
pub(crate) enum Value<'a> {
    Name(&'a str),
    Symbol(&'a str),
}

fn is_not_whitespace(c: char) -> bool {
    !c.is_whitespace()
}

fn name(input: &str) -> IResult<&str, Value<'_>> {
    map(all_consuming(take_while(is_not_whitespace)), Value::Name)(input)
}

fn symbol(input: &str) -> IResult<&str, Value<'_>> {
    map(
        delimited(tag("<"), take_until(">"), tag(">")),
        Value::Symbol,
    )(input)
}

pub(crate) fn parse(input: &str) -> IResult<&str, Value<'_>> {
    alt((symbol, name))(input)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() {
        assert_eq!(parse("sq0"), Ok(("", Value::Name("sq0"))));
        assert_eq!(parse("<sq0>"), Ok(("", Value::Symbol("sq0"))));

        assert_eq!(
            parse("sq 0"),
            Err(nom::Err::Error(nom::error::Error::new(
                " 0",
                nom::error::ErrorKind::Eof
            )))
        );
        assert_eq!(
            parse("sq\u{a0}0"),
            Err(nom::Err::Error(nom::error::Error::new(
                "\u{a0}0",
                nom::error::ErrorKind::Eof
            )))
        );
    }
}
