use std::{error, fmt, str::FromStr};

/// A VCF number describing the cardinality of a field.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Number {
    /// An explicit size.
    Count(usize),
    /// The number of alternate alleles.
    A,
    /// The number of reference and alternate alleles.
    R,
    /// The number of genotypes.
    G,
    /// The size is unknown.
    Unknown,
}

impl Default for Number {
    fn default() -> Self {
        Self::Count(1)
    }
}

impl fmt::Display for Number {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Count(n) => write!(f, "{n}"),
            Self::A => f.write_str("A"),
            Self::R => f.write_str("R"),
            Self::G => f.write_str("G"),
            Self::Unknown => f.write_str("."),
        }
    }
}

/// An error returned when a VCF number fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid => f.write_str("invalid input"),
        }
    }
}

impl FromStr for Number {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            "A" => Ok(Self::A),
            "R" => Ok(Self::R),
            "G" => Ok(Self::G),
            "." => Ok(Self::Unknown),
            _ => match s.parse() {
                Ok(n) => Ok(Self::Count(n)),
                Err(_) => Err(ParseError::Invalid),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(Number::default(), Number::Count(1));
    }

    #[test]
    fn test_fmt() {
        assert_eq!(Number::Count(1).to_string(), "1");
        assert_eq!(Number::A.to_string(), "A");
        assert_eq!(Number::R.to_string(), "R");
        assert_eq!(Number::G.to_string(), "G");
        assert_eq!(Number::Unknown.to_string(), ".");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("1".parse(), Ok(Number::Count(1)));
        assert_eq!("A".parse(), Ok(Number::A));
        assert_eq!("R".parse(), Ok(Number::R));
        assert_eq!("G".parse(), Ok(Number::G));
        assert_eq!(".".parse(), Ok(Number::Unknown));

        assert_eq!("".parse::<Number>(), Err(ParseError::Empty));
        assert_eq!("Noodles".parse::<Number>(), Err(ParseError::Invalid));
    }
}
