use std::{error, fmt, str::FromStr};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Number {
    Count(usize),
    A,
    R,
    G,
    Unknown,
}

impl Default for Number {
    fn default() -> Self {
        Self::Unknown
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid info number: expected {{<usize>, A, R, G, .}}, got {}",
            self.0
        )
    }
}

impl FromStr for Number {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "A" => Ok(Self::A),
            "R" => Ok(Self::R),
            "G" => Ok(Self::G),
            "." => Ok(Self::Unknown),
            _ => match s.parse() {
                Ok(n) => Ok(Self::Count(n)),
                Err(_) => Err(ParseError(s.into())),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(Number::default(), Number::Unknown);
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("1".parse::<Number>()?, Number::Count(1));
        assert_eq!("A".parse::<Number>()?, Number::A);
        assert_eq!("R".parse::<Number>()?, Number::R);
        assert_eq!("G".parse::<Number>()?, Number::G);
        assert_eq!(".".parse::<Number>()?, Number::Unknown);

        assert!("".parse::<Number>().is_err());
        assert!("Noodles".parse::<Number>().is_err());

        Ok(())
    }
}
