//! VCF header pedigree record.

pub mod key;

pub use self::key::Key;

use std::{convert::TryFrom, error, fmt};

use indexmap::IndexMap;

use super::{record, Record};

/// A VCF header pedigree record (`PEDIGREE`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Pedigree {
    id: String,
    fields: IndexMap<String, String>,
}

impl Pedigree {
    /// Creates a VCF header pedigree record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Pedigree;
    /// let pedigree = Pedigree::new(String::from("cid"), Default::default());
    /// ```
    pub fn new(id: String, fields: IndexMap<String, String>) -> Self {
        Self { id, fields }
    }

    /// Returns the ID of the sample record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Pedigree;
    /// let pedigree = Pedigree::new(String::from("cid"), Default::default());
    /// assert_eq!(pedigree.id(), "cid");
    /// ```
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Returns the extra fields in the record.
    ///
    /// This includes fields other than `ID`.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Pedigree;
    /// let pedigree = Pedigree::new(String::from("cid"), Default::default());
    /// assert!(pedigree.fields().is_empty());
    /// ```
    pub fn fields(&self) -> &IndexMap<String, String> {
        &self.fields
    }
}

impl fmt::Display for Pedigree {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(record::PREFIX)?;
        f.write_str(record::Key::Pedigree.as_ref())?;
        f.write_str("=<")?;

        write!(f, "ID={}", self.id())?;

        for (key, value) in &self.fields {
            write!(f, r#",{}={}"#, key, value)?;
        }

        f.write_str(">")?;

        Ok(())
    }
}

/// An error returned when a generic VCF header record fails to convert to a pedigree record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromRecordError {
    /// The record is invalid.
    InvalidRecord,
    /// A required field is missing.
    MissingField(Key),
}

impl error::Error for TryFromRecordError {}

impl fmt::Display for TryFromRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord => f.write_str("invalid record"),
            Self::MissingField(key) => write!(f, "missing field: {}", key),
        }
    }
}

impl TryFrom<Record> for Pedigree {
    type Error = TryFromRecordError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        match record.into() {
            (record::Key::Pedigree, record::Value::Struct(fields)) => parse_struct(fields),
            _ => Err(TryFromRecordError::InvalidRecord),
        }
    }
}

fn parse_struct(fields: Vec<(String, String)>) -> Result<Pedigree, TryFromRecordError> {
    let mut it = fields.into_iter();

    let id = it
        .next()
        .ok_or(TryFromRecordError::MissingField(Key::Id))
        .and_then(|(k, v)| match k.parse() {
            Ok(Key::Id) => Ok(v),
            _ => Err(TryFromRecordError::MissingField(Key::Id)),
        })?;

    let fields = it.collect();

    Ok(Pedigree::new(id, fields))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_record() -> Record {
        Record::new(
            record::Key::Pedigree,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("cid")),
                (String::from("Father"), String::from("fid")),
                (String::from("Mother"), String::from("mid")),
            ]),
        )
    }

    #[test]
    fn test_fmt() {
        let sample = Pedigree::new(String::from("cid"), IndexMap::new());
        assert_eq!(sample.to_string(), "##PEDIGREE=<ID=cid>");

        let mut fields = IndexMap::new();
        fields.insert(String::from("Father"), String::from("fid"));
        fields.insert(String::from("Mother"), String::from("mid"));
        let sample = Pedigree::new(String::from("cid"), fields);
        assert_eq!(
            sample.to_string(),
            "##PEDIGREE=<ID=cid,Father=fid,Mother=mid>"
        );
    }

    #[test]
    fn test_try_from_record_for_pedigree() {
        let record = build_record();
        let actual = Pedigree::try_from(record);

        let mut fields = IndexMap::new();
        fields.insert(String::from("Father"), String::from("fid"));
        fields.insert(String::from("Mother"), String::from("mid"));
        let expected = Pedigree::new(String::from("cid"), fields);

        assert_eq!(actual, Ok(expected));
    }
}
