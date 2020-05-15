use std::{error, fmt, str::FromStr};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Type {
    Integer,
    Float,
    Flag,
    Character,
    String,
}

impl AsRef<str> for Type {
    fn as_ref(&self) -> &str {
        match self {
            Self::Integer => "Integer",
            Self::Float => "Float",
            Self::Flag => "Flag",
            Self::Character => "Character",
            Self::String => "String",
        }
    }
}

impl Default for Type {
    fn default() -> Self {
        Self::String
    }
}

impl fmt::Display for Type {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid info type: expected {{Integer, Float, Flag, Character, String}}, got {}",
            self.0
        )
    }
}

impl FromStr for Type {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Integer" => Ok(Self::Integer),
            "Float" => Ok(Self::Float),
            "Flag" => Ok(Self::Flag),
            "Character" => Ok(Self::Character),
            "String" => Ok(Self::String),
            _ => Err(ParseError(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(Type::default(), Type::String);
    }

    #[test]
    fn test_fmt() {
        assert_eq!(Type::Integer.to_string(), "Integer");
        assert_eq!(Type::Float.to_string(), "Float");
        assert_eq!(Type::Flag.to_string(), "Flag");
        assert_eq!(Type::Character.to_string(), "Character");
        assert_eq!(Type::String.to_string(), "String");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("Integer".parse::<Type>()?, Type::Integer);
        assert_eq!("Float".parse::<Type>()?, Type::Float);
        assert_eq!("Flag".parse::<Type>()?, Type::Flag);
        assert_eq!("Character".parse::<Type>()?, Type::Character);
        assert_eq!("String".parse::<Type>()?, Type::String);

        assert!("".parse::<Type>().is_err());
        assert!("Noodles".parse::<Type>().is_err());

        Ok(())
    }
}
