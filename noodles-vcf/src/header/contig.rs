//! VCF header contig record and key.

use std::{error, fmt, num};

use indexmap::IndexMap;

use super::{record, Record};
use crate::record::chromosome;

const ID: &str = "ID";
const LENGTH: &str = "length";
const IDX: &str = "IDX";

/// A VCF header contig record (`contig`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Contig {
    id: String,
    len: Option<i32>,
    idx: Option<usize>,
    fields: IndexMap<String, String>,
}

#[allow(clippy::len_without_is_empty)]
impl Contig {
    /// Creates a VCF header contig record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Contig;
    /// let contig = Contig::new("sq0");
    /// ```
    pub fn new<S>(id: S) -> Self
    where
        S: Into<String>,
    {
        Self {
            id: id.into(),
            len: None,
            idx: None,
            fields: IndexMap::new(),
        }
    }

    /// Returns the ID of the contig.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Contig;
    /// let contig = Contig::new("sq0");
    /// assert_eq!(contig.id(), "sq0");
    /// ```
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Returns the length of the contig, if it is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Contig;
    /// let contig = Contig::new("sq0");
    /// assert_eq!(contig.len(), None);
    /// ```
    pub fn len(&self) -> Option<i32> {
        self.len
    }

    /// Returns a mutable reference to the length of the contig.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Contig;
    ///
    /// let mut contig = Contig::new("sq0");
    /// assert!(contig.len().is_none());
    ///
    /// *contig.len_mut() = Some(8);
    /// assert_eq!(contig.len(), Some(8));
    /// ```
    pub fn len_mut(&mut self) -> &mut Option<i32> {
        &mut self.len
    }

    /// Returns the index of the ID in the dictionary of strings.
    ///
    /// This is typically used in BCF.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Contig;
    /// let contig = Contig::new("sq0");
    /// assert!(contig.idx().is_none());
    /// ```
    pub fn idx(&self) -> Option<usize> {
        self.idx
    }

    /// Returns the value of the field with the given key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::{contig, record, Record, Contig};
    ///
    /// let record = Record::new(
    ///     record::Key::Contig,
    ///     record::Value::Struct(vec![
    ///         (String::from("ID"), String::from("sq0")),
    ///         (String::from("md5"), String::from("d7eba311421bbc9d3ada44709dd61534")),
    ///     ]),
    /// );
    /// let contig = Contig::try_from(record)?;
    ///
    /// assert_eq!(contig.get("md5"), Some("d7eba311421bbc9d3ada44709dd61534"));
    /// assert!(contig.get("species").is_none());
    ///
    /// # Ok::<(), contig::TryFromRecordError>(())
    /// ```
    pub fn get(&self, key: &str) -> Option<&str> {
        self.fields.get(key).map(|s| &**s)
    }
}

impl fmt::Display for Contig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(record::PREFIX)?;
        f.write_str(record::Key::Contig.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", ID, self.id)?;

        if let Some(len) = self.len {
            write!(f, ",{}={}", LENGTH, len)?;
        }

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

/// An error returned when a generic VCF header record fails to convert to a contig header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromRecordError {
    /// The record is invalid.
    InvalidRecord,
    /// The ID is invalid.
    InvalidId,
    /// The length is invalid.
    InvalidLength(num::ParseIntError),
    /// A required field is missing.
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
            Self::InvalidId => f.write_str("invalid ID"),
            Self::InvalidLength(e) => write!(f, "invalid length: {}", e),
            Self::InvalidIdx(e) => write!(f, "invalid index (`{}`): {}", IDX, e),
        }
    }
}

impl TryFrom<Record> for Contig {
    type Error = TryFromRecordError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        match record.into() {
            (record::Key::Contig, record::Value::Struct(fields)) => parse_struct(fields),
            _ => Err(TryFromRecordError::InvalidRecord),
        }
    }
}

fn parse_struct(fields: Vec<(String, String)>) -> Result<Contig, TryFromRecordError> {
    let mut id = None;
    let mut len = None;
    let mut idx = None;
    let mut other_fields = IndexMap::new();

    for (key, value) in fields {
        match key.as_ref() {
            ID => {
                if !chromosome::is_valid_name(&value) {
                    return Err(TryFromRecordError::InvalidId);
                }

                id = Some(value);
            }
            LENGTH => {
                len = value
                    .parse()
                    .map(Some)
                    .map_err(TryFromRecordError::InvalidLength)?;
            }
            IDX => {
                idx = value
                    .parse()
                    .map(Some)
                    .map_err(TryFromRecordError::InvalidIdx)?;
            }
            _ => {
                other_fields.insert(key, value);
            }
        }
    }

    Ok(Contig {
        id: id.ok_or(TryFromRecordError::MissingField(ID))?,
        len,
        idx,
        fields: other_fields,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_record() -> Record {
        Record::new(
            record::Key::Contig,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("sq0")),
                (String::from("length"), String::from("13")),
                (
                    String::from("md5"),
                    String::from("d7eba311421bbc9d3ada44709dd61534"),
                ),
            ]),
        )
    }

    #[test]
    fn test_fmt() -> Result<(), TryFromRecordError> {
        let record = build_record();
        let contig = Contig::try_from(record)?;

        let expected = r#"##contig=<ID=sq0,length=13,md5="d7eba311421bbc9d3ada44709dd61534">"#;
        assert_eq!(contig.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_contig() {
        let record = build_record();

        assert_eq!(
            Contig::try_from(record),
            Ok(Contig {
                id: String::from("sq0"),
                len: Some(13),
                idx: None,
                fields: [(
                    String::from("md5"),
                    String::from("d7eba311421bbc9d3ada44709dd61534")
                )]
                .into_iter()
                .collect(),
            })
        );
    }

    #[test]
    fn test_try_from_record_for_contig_with_extra_fields() {
        let record = Record::new(
            record::Key::Contig,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("sq0")),
                (String::from("length"), String::from("13")),
                (
                    String::from("md5"),
                    String::from("d7eba311421bbc9d3ada44709dd61534"),
                ),
                (String::from("IDX"), String::from("1")),
            ]),
        );

        assert_eq!(
            Contig::try_from(record),
            Ok(Contig {
                id: String::from("sq0"),
                len: Some(13),
                idx: Some(1),
                fields: [(
                    String::from("md5"),
                    String::from("d7eba311421bbc9d3ada44709dd61534")
                )]
                .into_iter()
                .collect(),
            })
        );
    }

    #[test]
    fn test_try_from_record_for_contig_with_an_invalid_record_key() {
        let record = Record::new(
            record::Key::FileFormat,
            record::Value::Struct(vec![(String::from("ID"), String::from("sq0"))]),
        );

        assert_eq!(
            Contig::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_contig_with_an_invalid_record_value() {
        let record = Record::new(
            record::Key::Contig,
            record::Value::String(String::from("VCF4.3")),
        );

        assert_eq!(
            Contig::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_contig_with_an_invalid_id() {
        let record = Record::new(
            record::Key::Contig,
            record::Value::Struct(vec![(String::from("ID"), String::from("sq 0"))]),
        );

        assert_eq!(Contig::try_from(record), Err(TryFromRecordError::InvalidId));
    }

    #[test]
    fn test_try_from_record_for_contig_with_an_invalid_length() {
        let record = Record::new(
            record::Key::Contig,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("sq0")),
                (String::from("length"), String::from("NA")),
            ]),
        );

        assert!(matches!(
            Contig::try_from(record),
            Err(TryFromRecordError::InvalidLength(_))
        ));
    }

    #[test]
    fn test() {
        let record = Record::new(
            record::Key::Contig,
            record::Value::Struct(vec![
                (String::from("ID"), String::from("sq0")),
                (String::from("length"), String::from("13")),
                (String::from("IDX"), String::from("ndls")),
            ]),
        );

        assert!(matches!(
            Contig::try_from(record),
            Err(TryFromRecordError::InvalidIdx(_))
        ));
    }
}
