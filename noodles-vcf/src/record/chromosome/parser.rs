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
    Reference(&'a str),
}

fn name(input: &str) -> IResult<&str, Value> {
    map(rest, |s| Value::Name(s))(input)
}

fn reference(input: &str) -> IResult<&str, Value> {
    map(delimited(tag("<"), take_until(">"), tag(">")), |s| {
        Value::Reference(s)
    })(input)
}

pub(crate) fn parse(input: &str) -> IResult<&str, Value> {
    alt((reference, name))(input)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() {
        assert_eq!(parse("sq0"), Ok(("", Value::Name("sq0"))));
        assert_eq!(parse("<sq0>"), Ok(("", Value::Reference("sq0"))));
    }
}
