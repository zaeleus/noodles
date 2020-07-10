pub mod ty;

pub use self::ty::Type;

use std::{error, fmt, str::FromStr};

const DELIMITER: char = ':';

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StructuralVariant {
    ty: Type,
    subtypes: Vec<String>,
}

impl StructuralVariant {
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

impl fmt::Display for StructuralVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.ty.as_ref())?;

        for subtype in self.subtypes() {
            write!(f, "{}{}", DELIMITER, subtype)?;
        }

        Ok(())
    }
}

impl From<Type> for StructuralVariant {
    fn from(ty: Type) -> Self {
        Self::new(ty, Vec::new())
    }
}

/// An error returned when a raw VCF record a structural variant of an alternate base allele symbol
/// fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The type is missing.
    MissingType,
    /// The type is invalid.
    InvalidType(ty::ParseError),
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::MissingType => f.write_str("missing type"),
            Self::InvalidType(e) => write!(f, "invalid type: {}", e),
        }
    }
}

impl error::Error for ParseError {}

impl FromStr for StructuralVariant {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

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
        let symbol = StructuralVariant::new(Type::Deletion, Vec::new());
        assert_eq!(symbol.to_string(), "DEL");

        let symbol = StructuralVariant::new(Type::Deletion, vec![String::from("TANDEM")]);
        assert_eq!(symbol.to_string(), "DEL:TANDEM");

        let symbol = StructuralVariant::new(
            Type::Deletion,
            vec![String::from("ME"), String::from("ALU")],
        );
        assert_eq!(symbol.to_string(), "DEL:ME:ALU");
    }

    #[test]
    fn test_from_type_for_symbol() {
        assert_eq!(
            StructuralVariant::from(Type::Deletion),
            StructuralVariant::new(Type::Deletion, Vec::new())
        );
    }

    #[test]
    fn test_from_str() {
        assert_eq!(
            "DEL".parse(),
            Ok(StructuralVariant::new(Type::Deletion, Vec::new()))
        );

        assert_eq!(
            "DEL:TANDEM".parse(),
            Ok(StructuralVariant::new(
                Type::Deletion,
                vec![String::from("TANDEM")]
            ))
        );

        assert_eq!(
            "DEL:ME:ALU".parse(),
            Ok(StructuralVariant::new(
                Type::Deletion,
                vec![String::from("ME"), String::from("ALU")]
            ))
        );

        assert_eq!("".parse::<StructuralVariant>(), Err(ParseError::Empty));
        assert!(matches!(
            "NDL".parse::<StructuralVariant>(),
            Err(ParseError::InvalidType(_))
        ));
    }
}
