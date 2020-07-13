//! VCF record alternate bases allele structural variant symbol and type.

pub mod ty;

pub use self::ty::Type;

use std::{error, fmt, str::FromStr};

const DELIMITER: char = ':';

/// A VCF alternate bases allele structural variant symbol.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StructuralVariant {
    ty: Type,
    subtypes: Vec<String>,
}

impl StructuralVariant {
    /// Creates a allele structural variant symbol.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::alternate_bases::allele::symbol::{
    ///     structural_variant::Type,
    ///     StructuralVariant,
    /// };
    ///
    /// let structural_variant = StructuralVariant::new(Type::Deletion, Vec::new());
    /// ```
    pub fn new(ty: Type, subtypes: Vec<String>) -> Self {
        Self { ty, subtypes }
    }

    /// Returns the structural variant type.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::alternate_bases::allele::symbol::{
    ///     structural_variant::Type,
    ///     StructuralVariant,
    /// };
    ///
    /// let structural_variant = StructuralVariant::new(Type::Deletion, Vec::new());
    /// assert_eq!(structural_variant.ty(), Type::Deletion);
    /// ```
    pub fn ty(&self) -> Type {
        self.ty
    }

    /// Returns the structural variant subtypes.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::alternate_bases::allele::symbol::{
    ///     structural_variant::Type,
    ///     StructuralVariant,
    /// };
    ///
    /// let structural_variant = StructuralVariant::new(
    ///     Type::Deletion,
    ///     vec![String::from("ME")],
    /// );
    ///
    /// assert_eq!(structural_variant.subtypes(), [String::from("ME")]);
    /// ```
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
