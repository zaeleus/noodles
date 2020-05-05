use std::{error, fmt, str::FromStr};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Key {
    Id,
    Number,
    Type,
    Description,
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid format key: expected {{ID, Number, Type, Description}}, got {}",
            self.0
        )
    }
}

impl FromStr for Key {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "ID" => Ok(Self::Id),
            "Number" => Ok(Self::Number),
            "Type" => Ok(Self::Type),
            "Description" => Ok(Self::Description),
            _ => Err(ParseError(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("ID".parse::<Key>()?, Key::Id);
        assert_eq!("Number".parse::<Key>()?, Key::Number);
        assert_eq!("Type".parse::<Key>()?, Key::Type);
        assert_eq!("Description".parse::<Key>()?, Key::Description);

        assert!("".parse::<Key>().is_err());
        assert!("Noodles".parse::<Key>().is_err());

        Ok(())
    }
}
