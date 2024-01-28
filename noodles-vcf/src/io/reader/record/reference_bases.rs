use std::{error, fmt};

/// An error when raw VCF record reference bases fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
        }
    }
}

pub(super) fn parse_reference_bases(
    s: &str,
    reference_bases: &mut String,
) -> Result<(), ParseError> {
    if s.is_empty() {
        return Err(ParseError::Empty);
    }

    reference_bases.clear();
    reference_bases.push_str(s);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_reference_bases() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = String::new();

        parse_reference_bases("ATCGN", &mut buf)?;
        assert_eq!(buf, "ATCGN");

        parse_reference_bases("atcgn", &mut buf)?;
        assert_eq!(buf, "atcgn");

        parse_reference_bases("AtCgN", &mut buf)?;
        assert_eq!(buf, "AtCgN");

        assert_eq!(parse_reference_bases("", &mut buf), Err(ParseError::Empty));

        Ok(())
    }
}
