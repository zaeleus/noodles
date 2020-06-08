use nom::{
    branch::alt,
    bytes::complete::{tag, take_until},
    combinator::{map, rest},
    sequence::delimited,
    IResult,
};

#[derive(Debug, Eq, PartialEq)]
pub(crate) enum Value<'a> {
    Name(&'a str),
    Symbol(&'a str),
}

fn name(input: &str) -> IResult<&str, Value<'_>> {
    map(rest, Value::Name)(input)
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
    }
}
