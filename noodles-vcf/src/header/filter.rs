//! VCF header filter record and key.

use std::{error, fmt, num};

use indexmap::IndexMap;

use super::{record, Record};
use crate::record::Filters;

const ID: &str = "ID";
const DESCRIPTION: &str = "Description";
const IDX: &str = "IDX";

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
    /// assert_eq!(filter, Filter::new("PASS", "All filters passed"));
    /// ```
    pub fn pass() -> Self {
        Self::new(Filters::Pass.to_string(), "All filters passed")
    }

    pub(crate) fn try_from_fields(
        id: String,
        fields: IndexMap<String, String>,
    ) -> Result<Self, TryFromRecordError> {
        parse_struct(id, fields)
    }

    /// Creates a VCF header filter record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Filter;
    /// let filter = Filter::new("q10", "Quality below 10");
    /// ```
    pub fn new<S, T>(id: S, description: T) -> Self
    where
        S: Into<String>,
        T: Into<String>,
    {
        Self {
            id: id.into(),
            description: description.into(),
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
    /// let filter = Filter::new("q10", "Quality below 10");
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
    /// let filter = Filter::new("q10", "Quality below 10");
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
    /// let filter = Filter::new("q10", "Quality below 10");
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
    /// let filter = Filter::new("q10", "Quality below 10");
    /// assert!(filter.fields().is_empty());
    /// ```
    pub fn fields(&self) -> &IndexMap<String, String> {
        &self.fields
    }
}

impl fmt::Display for Filter {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(record::PREFIX)?;
        f.write_str(record::key::FILTER.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", ID, self.id)?;

        write!(f, ",{}=", DESCRIPTION)?;
        super::fmt::write_escaped_string(f, self.description())?;

        for (key, value) in &self.fields {
            write!(f, ",{}=", key)?;
            super::fmt::write_escaped_string(f, value)?;
        }

        if let Some(idx) = self.idx() {
            write!(f, ",{}={}", IDX, idx)?;
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
    MissingField(&'static str),
    /// The index (`IDX`) is invalid.
    InvalidIdx(num::ParseIntError),
}

impl error::Error for TryFromRecordError {}

impl fmt::Display for TryFromRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord => f.write_str("invalid record"),
            Self::MissingField(key) => write!(f, "missing field: {}", key),
            Self::InvalidIdx(e) => write!(f, "invalid index (`{}`): {}", IDX, e),
        }
    }
}

impl TryFrom<Record> for Filter {
    type Error = TryFromRecordError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        match record.into() {
            (record::key::FILTER, record::Value::Struct(id, fields)) => parse_struct(id, fields),
            _ => Err(TryFromRecordError::InvalidRecord),
        }
    }
}

fn parse_struct(
    id: String,
    fields: IndexMap<String, String>,
) -> Result<Filter, TryFromRecordError> {
    let mut it = fields.into_iter();

    let description = it
        .next()
        .ok_or(TryFromRecordError::MissingField(DESCRIPTION))
        .and_then(|(k, v)| match k.as_ref() {
            DESCRIPTION => Ok(v),
            _ => Err(TryFromRecordError::MissingField(DESCRIPTION)),
        })?;

    let mut idx = None;
    let mut fields = IndexMap::new();

    for (key, value) in it {
        match key.as_ref() {
            IDX => {
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
            record::key::FILTER,
            record::Value::Struct(
                String::from("q10"),
                [(
                    String::from("Description"),
                    String::from("Quality below 10"),
                )]
                .into_iter()
                .collect(),
            ),
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
            Ok(Filter::new("q10", "Quality below 10"))
        );
    }

    #[test]
    fn test_try_from_record_for_filter_with_extra_fields() {
        let record = Record::new(
            record::key::FILTER,
            record::Value::Struct(
                String::from("q10"),
                [
                    (
                        String::from("Description"),
                        String::from("Quality below 10"),
                    ),
                    (String::from("Source"), String::from("noodles")),
                    (String::from("IDX"), String::from("1")),
                ]
                .into_iter()
                .collect(),
            ),
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
            record::key::FILE_FORMAT,
            record::Value::Struct(
                String::from("q10"),
                [(
                    String::from("Description"),
                    String::from("Quality below 10"),
                )]
                .into_iter()
                .collect(),
            ),
        );

        assert_eq!(
            Filter::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_filter_with_an_invalid_record_value() {
        let record = Record::new(
            record::key::FILTER,
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
            record::key::FILTER,
            record::Value::Struct(String::from("q10"), Default::default()),
        );

        assert!(matches!(
            Filter::try_from(record),
            Err(TryFromRecordError::MissingField(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_filter_with_an_invalid_idx() {
        let record = Record::new(
            record::key::FILTER,
            record::Value::Struct(
                String::from("q10"),
                [
                    (
                        String::from("Description"),
                        String::from("Quality below 10"),
                    ),
                    (String::from("IDX"), String::from("ndls")),
                ]
                .into_iter()
                .collect(),
            ),
        );

        assert!(matches!(
            Filter::try_from(record),
            Err(TryFromRecordError::InvalidIdx(_))
        ));
    }
}
