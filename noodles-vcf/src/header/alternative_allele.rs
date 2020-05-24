mod key;

use std::{convert::TryFrom, error, fmt};

use crate::record::alternate_bases::allele::{symbol, Symbol};

use super::record;

use self::key::Key;

#[derive(Clone, Debug)]
pub struct AlternativeAllele {
    id: Symbol,
    description: String,
}

impl AlternativeAllele {
    pub fn new(id: Symbol, description: String) -> Self {
        Self { id, description }
    }

    pub fn id(&self) -> &Symbol {
        &self.id
    }

    pub fn description(&self) -> &str {
        &self.description
    }
}

impl fmt::Display for AlternativeAllele {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("##")?;
        f.write_str(record::Kind::AlternativeAllele.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", Key::Id, self.id)?;
        write!(f, r#",{}="{}""#, Key::Description, self.description)?;

        f.write_str(">")?;

        Ok(())
    }
}

#[derive(Debug)]
pub enum ParseError {
    MissingField(Key),
    InvalidId(symbol::ParseError),
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid alternative allele header: ")?;

        match self {
            ParseError::MissingField(key) => write!(f, "missing {} field", key),
            ParseError::InvalidId(e) => write!(f, "{}", e),
        }
    }
}

impl error::Error for ParseError {}

impl TryFrom<&[(String, String)]> for AlternativeAllele {
    type Error = ParseError;

    fn try_from(fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut it = fields.iter();

        let id = it
            .next()
            .ok_or_else(|| ParseError::MissingField(Key::Id))
            .and_then(|(k, v)| match k.parse() {
                Ok(Key::Id) => v.parse().map_err(ParseError::InvalidId),
                _ => Err(ParseError::MissingField(Key::Id)),
            })?;

        let description = it
            .next()
            .ok_or_else(|| ParseError::MissingField(Key::Description))
            .and_then(|(k, v)| match k.parse() {
                Ok(Key::Description) => Ok(v.into()),
                _ => Err(ParseError::MissingField(Key::Description)),
            })?;

        Ok(Self { id, description })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_fields() -> Vec<(String, String)> {
        vec![
            (String::from("ID"), String::from("DEL")),
            (String::from("Description"), String::from("Deletion")),
        ]
    }

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let fields = build_fields();
        let alternative_allele = AlternativeAllele::try_from(&fields[..])?;

        let expected = r#"##ALT=<ID=DEL,Description="Deletion">"#;

        assert_eq!(alternative_allele.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_filter() -> Result<(), ParseError> {
        let fields = build_fields();
        let filter = AlternativeAllele::try_from(&fields[..])?;

        assert_eq!(
            filter.id(),
            &Symbol::StructuralVariant(symbol::StructuralVariant::from(
                symbol::structural_variant::Type::Deletion
            ))
        );
        assert_eq!(filter.description(), "Deletion");

        Ok(())
    }
}
