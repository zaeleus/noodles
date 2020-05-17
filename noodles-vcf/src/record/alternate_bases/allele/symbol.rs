mod ty;

pub use self::ty::Type;

use std::{error, fmt, str::FromStr};

const DELIMITER: char = ':';

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Symbol {
    ty: Type,
    subtypes: Vec<String>,
}

impl Symbol {
    pub fn new(ty: Type, subtypes: Vec<String>) -> Self {
        Self { ty, subtypes }
    }

    pub fn ty(&self) -> Type {
        self.ty
    }

    pub fn subtypes(&self) -> &[String] {
        &self.subtypes
    }
}

impl fmt::Display for Symbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.ty.as_ref())?;

        for subtype in self.subtypes() {
            write!(f, "{}{}", DELIMITER, subtype)?;
        }

        Ok(())
    }
}

impl From<Type> for Symbol {
    fn from(ty: Type) -> Self {
        Self::new(ty, Vec::new())
    }
}

#[derive(Debug)]
pub enum ParseError {
    MissingType,
    InvalidType(ty::ParseError),
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid alternate bases allele symbol")?;

        match self {
            Self::MissingType => f.write_str("missing type"),
            Self::InvalidType(e) => write!(f, "{}", e),
        }
    }
}

impl error::Error for ParseError {}

impl FromStr for Symbol {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut components = s.split(DELIMITER);

        let ty = components
            .next()
            .ok_or_else(|| ParseError::MissingType)
            .and_then(|s| s.parse().map_err(ParseError::InvalidType))?;

        let subtypes = components.map(|s| s.into()).collect();

        Ok(Self::new(ty, subtypes))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let symbol = Symbol::new(Type::Deletion, Vec::new());
        assert_eq!(symbol.to_string(), "DEL");

        let symbol = Symbol::new(Type::Deletion, vec![String::from("TANDEM")]);
        assert_eq!(symbol.to_string(), "DEL:TANDEM");

        let symbol = Symbol::new(
            Type::Deletion,
            vec![String::from("ME"), String::from("ALU")],
        );
        assert_eq!(symbol.to_string(), "DEL:ME:ALU");
    }

    #[test]
    fn test_from_type_for_symbol() {
        let symbol = Symbol::from(Type::Deletion);
        assert_eq!(symbol.ty(), Type::Deletion);
        assert!(symbol.subtypes().is_empty());
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let actual: Symbol = "DEL".parse()?;
        assert_eq!(actual.ty(), Type::Deletion);
        assert!(actual.subtypes().is_empty());

        let actual: Symbol = "DEL:TANDEM".parse()?;
        assert_eq!(actual.ty(), Type::Deletion);
        assert_eq!(actual.subtypes(), &[String::from("TANDEM")]);

        let actual: Symbol = "DEL:ME:ALU".parse()?;
        assert_eq!(actual.ty(), Type::Deletion);
        assert_eq!(
            actual.subtypes(),
            &[String::from("ME"), String::from("ALU")]
        );

        assert!("".parse::<Symbol>().is_err());

        Ok(())
    }
}
