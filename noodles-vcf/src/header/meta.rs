//! VCF header meta record.

use std::{error, fmt};

use super::{record, Number, Record};

const ID: &str = "ID";
const TYPE: &str = "Type";
const NUMBER: &str = "Number";
const VALUES: &str = "Values";

/// A VCF header meta record (`META`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Meta {
    id: String,
    values: Vec<String>,
}

impl Meta {
    /// Creates a VCF header meta record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Meta;
    ///
    /// let meta = Meta::new(
    ///     String::from("Assay"),
    ///     vec![String::from("WholeGenome"), String::from("Exome")],
    /// );
    /// ```
    pub fn new(id: String, values: Vec<String>) -> Self {
        Self { id, values }
    }

    /// Returns the ID of the meta record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Meta;
    ///
    /// let meta = Meta::new(
    ///     String::from("Assay"),
    ///     vec![String::from("WholeGenome"), String::from("Exome")],
    /// );
    ///
    /// assert_eq!(meta.id(), "Assay");
    /// ```
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Returns the values of the meta record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Meta;
    ///
    /// let meta = Meta::new(
    ///     String::from("Assay"),
    ///     vec![String::from("WholeGenome"), String::from("Exome")],
    /// );
    ///
    /// assert_eq!(meta.values(), ["WholeGenome", "Exome"]);
    /// ```
    pub fn values(&self) -> &[String] {
        &self.values
    }
}

impl fmt::Display for Meta {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(record::PREFIX)?;
        f.write_str(record::Key::Meta.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", ID, self.id)?;
        write!(f, ",{}=String", TYPE)?;
        write!(f, ",{}={}", NUMBER, Number::Unknown)?;

        write!(f, ",{}=[", VALUES)?;

        for (i, value) in self.values().iter().enumerate() {
            if i > 0 {
                f.write_str(", ")?;
            }

            f.write_str(value)?;
        }

        f.write_str("]")?;

        f.write_str(">")?;

        Ok(())
    }
}

/// An error returned when a generic VCF header record fails to convert to a meta header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromRecordError {
    /// The record is invalid.
    InvalidRecord,
    /// A required field is missing.
    MissingField(&'static str),
    /// The ID is invalid.
    InvalidId,
    /// The number is invalid.
    InvalidType,
    /// The type is invalid.
    InvalidNumber,
}

impl error::Error for TryFromRecordError {}

impl fmt::Display for TryFromRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord => f.write_str("invalid record"),
            Self::MissingField(key) => write!(f, "missing field: {}", key),
            Self::InvalidId => f.write_str("invalid ID"),
            Self::InvalidNumber => f.write_str("invalid number"),
            Self::InvalidType => f.write_str("invalid type"),
        }
    }
}

impl TryFrom<Record> for Meta {
    type Error = TryFromRecordError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        match record.into() {
            (record::Key::Meta, record::Value::Struct(fields)) => parse_struct(fields),
            _ => Err(TryFromRecordError::InvalidRecord),
        }
    }
}

fn parse_struct(fields: Vec<(String, String)>) -> Result<Meta, TryFromRecordError> {
    let mut it = fields.into_iter();

    let id = it
        .next()
        .ok_or(TryFromRecordError::MissingField(ID))
        .and_then(|(k, v)| match k.as_ref() {
            ID => Ok(v),
            _ => Err(TryFromRecordError::MissingField(ID)),
        })?;

    let ty = it
        .next()
        .ok_or(TryFromRecordError::MissingField(TYPE))
        .and_then(|(k, v)| match k.as_ref() {
            TYPE => Ok(v),
            _ => Err(TryFromRecordError::MissingField(TYPE)),
        })?;

    if ty != "String" {
        return Err(TryFromRecordError::InvalidType);
    }

    let number = it
        .next()
        .ok_or(TryFromRecordError::MissingField(NUMBER))
        .and_then(|(k, v)| match k.as_ref() {
            NUMBER => Ok(v),
            _ => Err(TryFromRecordError::MissingField(NUMBER)),
        })?;

    if number != "." {
        return Err(TryFromRecordError::InvalidNumber);
    }

    let values = it
        .next()
        .ok_or(TryFromRecordError::MissingField(VALUES))
        .and_then(|(k, v)| match k.as_ref() {
            VALUES => Ok(v.split(',').map(|s| s.trim().into()).collect()),
            _ => Err(TryFromRecordError::MissingField(VALUES)),
        })?;

    Ok(Meta::new(id, values))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_record() -> Record {
        Record::new(
            record::Key::Meta,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("Assay")),
                (String::from("Type"), String::from("String")),
                (String::from("Number"), String::from(".")),
                (String::from("Values"), String::from("WholeGenome, Exome")),
            ]),
        )
    }

    #[test]
    fn test_fmt() -> Result<(), TryFromRecordError> {
        let record = build_record();
        let meta = Meta::try_from(record)?;

        let expected = "##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>";
        assert_eq!(meta.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_meta() {
        let record = build_record();
        let actual = Meta::try_from(record);

        let expected = Meta::new(
            String::from("Assay"),
            vec![String::from("WholeGenome"), String::from("Exome")],
        );

        assert_eq!(actual, Ok(expected));
    }
}
