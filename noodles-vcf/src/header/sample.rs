//! VCF header sample record.

use std::{error, fmt};

use indexmap::IndexMap;

use super::{record, Record};

const ID: &str = "ID";

/// A VCF header sample record (`SAMPLE`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Sample {
    id: String,
    fields: IndexMap<String, String>,
}

impl Sample {
    /// Creates a VCF header sample record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Sample;
    /// let sample = Sample::new(String::from("sample0"), Default::default());
    /// ```
    pub fn new(id: String, fields: IndexMap<String, String>) -> Self {
        Self { id, fields }
    }

    /// Returns the ID of the sample record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Sample;
    /// let sample = Sample::new(String::from("sample0"), Default::default());
    /// assert_eq!(sample.id(), "sample0");
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
    /// use noodles_vcf::header::Sample;
    /// let sample = Sample::new(String::from("sample0"), Default::default());
    /// assert!(sample.fields().is_empty());
    /// ```
    pub fn fields(&self) -> &IndexMap<String, String> {
        &self.fields
    }
}

impl fmt::Display for Sample {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(record::PREFIX)?;
        f.write_str(record::key::SAMPLE.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", ID, self.id())?;

        for (key, value) in &self.fields {
            write!(f, r#",{}={}"#, key, value)?;
        }

        f.write_str(">")?;

        Ok(())
    }
}

/// An error returned when a generic VCF header record fails to convert to a sample header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromRecordError {
    /// The record is invalid.
    InvalidRecord,
    /// A required field is missing.
    MissingField(&'static str),
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

impl TryFrom<Record> for Sample {
    type Error = TryFromRecordError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        match record.into() {
            (record::key::SAMPLE, record::Value::Struct(fields)) => parse_struct(fields),
            _ => Err(TryFromRecordError::InvalidRecord),
        }
    }
}

fn parse_struct(fields: Vec<(String, String)>) -> Result<Sample, TryFromRecordError> {
    let mut it = fields.into_iter();

    let id = it
        .next()
        .ok_or(TryFromRecordError::MissingField(ID))
        .and_then(|(k, v)| match k.as_ref() {
            ID => Ok(v),
            _ => Err(TryFromRecordError::MissingField(ID)),
        })?;

    let fields = it.collect();

    Ok(Sample::new(id, fields))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_record() -> Record {
        Record::new(
            record::key::SAMPLE,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("sample0")),
                (String::from("Assay"), String::from("WholeGenome")),
            ]),
        )
    }

    #[test]
    fn test_fmt() {
        let sample = Sample::new(String::from("sample0"), IndexMap::new());
        assert_eq!(sample.to_string(), "##SAMPLE=<ID=sample0>");

        let mut fields = IndexMap::new();
        fields.insert(String::from("Assay"), String::from("WholeGenome"));
        let sample = Sample::new(String::from("sample0"), fields);
        assert_eq!(
            sample.to_string(),
            "##SAMPLE=<ID=sample0,Assay=WholeGenome>"
        );
    }

    #[test]
    fn test_try_from_record_for_sample() {
        let record = build_record();
        let actual = Sample::try_from(record);

        let mut fields = IndexMap::new();
        fields.insert(String::from("Assay"), String::from("WholeGenome"));
        let expected = Sample::new(String::from("sample0"), fields);

        assert_eq!(actual, Ok(expected));
    }
}
