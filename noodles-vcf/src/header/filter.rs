mod key;

use std::{convert::TryFrom, error, fmt};

use super::record;

use self::key::Key;

/// A VCF header filter record (`FILTER`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Filter {
    id: String,
    description: String,
}

impl Filter {
    /// Creates a VCF header filter record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Filter;
    /// let filter = Filter::new(String::from("q10"), String::from("Quality below 10"));
    /// ```
    pub fn new(id: String, description: String) -> Self {
        Self { id, description }
    }

    /// Returns the ID of the filter.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Filter;
    /// let filter = Filter::new(String::from("q10"), String::from("Quality below 10"));
    /// assert_eq!(filter.id(), "q10");
    /// ```
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Returns the description of the filter.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Filter;
    /// let filter = Filter::new(String::from("q10"), String::from("Quality below 10"));
    /// assert_eq!(filter.description(), "Quality below 10");
    /// ```
    pub fn description(&self) -> &str {
        &self.description
    }
}

impl fmt::Display for Filter {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("##")?;
        f.write_str(record::Kind::Filter.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", Key::Id, self.id)?;
        write!(f, r#",{}="{}""#, Key::Description, self.description)?;

        f.write_str(">")?;

        Ok(())
    }
}

/// An error returned when a raw VCF header filter record fails to parse.
#[derive(Debug)]
pub enum ParseError {
    /// A field is missing.
    MissingField(Key),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid filter header: ")?;

        match self {
            ParseError::MissingField(key) => write!(f, "missing {} field", key),
        }
    }
}

impl TryFrom<&[(String, String)]> for Filter {
    type Error = ParseError;

    fn try_from(fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut it = fields.iter();

        let id = it
            .next()
            .ok_or_else(|| ParseError::MissingField(Key::Id))
            .and_then(|(k, v)| match k.parse() {
                Ok(Key::Id) => Ok(v.into()),
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
            (String::from("ID"), String::from("q10")),
            (
                String::from("Description"),
                String::from("Quality below 10"),
            ),
        ]
    }

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let fields = build_fields();
        let filter = Filter::try_from(&fields[..])?;

        let expected = r#"##FILTER=<ID=q10,Description="Quality below 10">"#;

        assert_eq!(filter.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_filter() -> Result<(), ParseError> {
        let fields = build_fields();
        let filter = Filter::try_from(&fields[..])?;

        assert_eq!(filter.id(), "q10");
        assert_eq!(filter.description(), "Quality below 10");

        Ok(())
    }
}
