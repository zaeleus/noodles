//! VCF header filter record and key.

pub mod key;

pub use self::key::Key;

use std::{convert::TryFrom, error, fmt, num};

use indexmap::IndexMap;

use crate::record::Filters;

use super::{record, Record};

/// A VCF header filter record (`FILTER`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Filter {
    id: String,
    description: String,
    idx: Option<usize>,
    fields: IndexMap<String, String>,
}

impl Filter {
    /// Creates a default filter record for PASS.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Filter;
    /// let filter = Filter::pass();
    /// assert_eq!(filter, Filter::new(String::from("PASS"), String::from("All filters passed")));
    /// ```
    pub fn pass() -> Self {
        Self::new(
            Filters::Pass.to_string(),
            String::from("All filters passed"),
        )
    }

    /// Creates a VCF header filter record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Filter;
    /// let filter = Filter::new(String::from("q10"), String::from("Quality below 10"));
    /// ```
    pub fn new(id: String, description: String) -> Self {
        Self {
            id,
            description,
            idx: None,
            fields: IndexMap::new(),
        }
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

    /// Returns the index of the ID in the dictionary of strings.
    ///
    /// This is typically used in BCF.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Filter;
    /// let filter = Filter::new(String::from("q10"), String::from("Quality below 10"));
    /// assert!(filter.idx().is_none());
    /// ```
    pub fn idx(&self) -> Option<usize> {
        self.idx
    }

    /// Returns the extra fields in the record.
    ///
    /// This includes fields other than `ID` and `Description`.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Filter;
    /// let filter = Filter::new(String::from("q10"), String::from("Quality below 10"));
    /// assert!(filter.fields().is_empty());
    /// ```
    pub fn fields(&self) -> &IndexMap<String, String> {
        &self.fields
    }
}

impl fmt::Display for Filter {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(record::PREFIX)?;
        f.write_str(record::Key::Filter.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", Key::Id, self.id)?;

        write!(f, ",{}=", Key::Description)?;
        super::fmt::write_escaped_string(f, self.description())?;

        if let Some(idx) = self.idx() {
            write!(f, ",{}={}", Key::Idx, idx)?;
        }

        for (key, value) in &self.fields {
            write!(f, ",{}=", key)?;
            super::fmt::write_escaped_string(f, value)?;
        }

        f.write_str(">")?;

        Ok(())
    }
}

/// An error returned when a generic VCF header record fails to convert to a filter header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromRecordError {
    /// The record is invalid.
    InvalidRecord,
    /// A field is missing.
    MissingField(Key),
    /// The index (`IDX`) is invalid.
    InvalidIdx(num::ParseIntError),
}

impl error::Error for TryFromRecordError {}

impl fmt::Display for TryFromRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord => f.write_str("invalid record"),
            Self::MissingField(key) => write!(f, "missing field: {}", key),
            Self::InvalidIdx(e) => write!(f, "invalid index (`{}`): {}", Key::Idx, e),
        }
    }
}

impl TryFrom<Record> for Filter {
    type Error = TryFromRecordError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        match record.into() {
            (record::Key::Filter, record::Value::Struct(fields)) => parse_struct(fields),
            _ => Err(TryFromRecordError::InvalidRecord),
        }
    }
}

fn parse_struct(fields: Vec<(String, String)>) -> Result<Filter, TryFromRecordError> {
    let mut it = fields.into_iter();

    let id = it
        .next()
        .ok_or(TryFromRecordError::MissingField(Key::Id))
        .and_then(|(k, v)| match k.parse() {
            Ok(Key::Id) => Ok(v),
            _ => Err(TryFromRecordError::MissingField(Key::Id)),
        })?;

    let description = it
        .next()
        .ok_or(TryFromRecordError::MissingField(Key::Description))
        .and_then(|(k, v)| match k.parse() {
            Ok(Key::Description) => Ok(v),
            _ => Err(TryFromRecordError::MissingField(Key::Description)),
        })?;

    let mut idx = None;
    let mut fields = IndexMap::new();

    for (key, value) in it {
        match key.parse() {
            Ok(Key::Idx) => {
                idx = value
                    .parse()
                    .map(Some)
                    .map_err(TryFromRecordError::InvalidIdx)?;
            }
            _ => {
                fields.insert(key, value);
            }
        }
    }

    Ok(Filter {
        id,
        description,
        idx,
        fields,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_record() -> Record {
        Record::new(
            record::Key::Filter,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("q10")),
                (
                    String::from("Description"),
                    String::from("Quality below 10"),
                ),
            ]),
        )
    }

    #[test]
    fn test_fmt() -> Result<(), TryFromRecordError> {
        let record = build_record();
        let filter = Filter::try_from(record)?;

        let expected = r#"##FILTER=<ID=q10,Description="Quality below 10">"#;
        assert_eq!(filter.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_filter() {
        let record = build_record();

        assert_eq!(
            Filter::try_from(record),
            Ok(Filter::new(
                String::from("q10"),
                String::from("Quality below 10"),
            ))
        );
    }

    #[test]
    fn test_try_from_record_for_filter_with_extra_fields() {
        let record = Record::new(
            record::Key::Filter,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("q10")),
                (
                    String::from("Description"),
                    String::from("Quality below 10"),
                ),
                (String::from("Source"), String::from("noodles")),
                (String::from("IDX"), String::from("1")),
            ]),
        );

        assert_eq!(
            Filter::try_from(record),
            Ok(Filter {
                id: String::from("q10"),
                description: String::from("Quality below 10"),
                idx: Some(1),
                fields: vec![(String::from("Source"), String::from("noodles"))]
                    .into_iter()
                    .collect()
            })
        );
    }

    #[test]
    fn test_try_from_record_for_filter_with_an_invalid_record_key() {
        let record = Record::new(
            record::Key::FileFormat,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("q10")),
                (
                    String::from("Description"),
                    String::from("Quality below 10"),
                ),
            ]),
        );

        assert_eq!(
            Filter::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_filter_with_an_invalid_record_value() {
        let record = Record::new(
            record::Key::Filter,
            record::Value::String(String::from("VCFv4.3")),
        );

        assert_eq!(
            Filter::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_filter_with_a_missing_field() {
        let record = Record::new(
            record::Key::Filter,
            record::Value::Struct(vec![(String::from("ID"), String::from("q10"))]),
        );

        assert!(matches!(
            Filter::try_from(record),
            Err(TryFromRecordError::MissingField(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_filter_with_an_invalid_idx() {
        let record = Record::new(
            record::Key::Filter,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("q10")),
                (
                    String::from("Description"),
                    String::from("Quality below 10"),
                ),
                (String::from("IDX"), String::from("ndls")),
            ]),
        );

        assert!(matches!(
            Filter::try_from(record),
            Err(TryFromRecordError::InvalidIdx(_))
        ));
    }
}
