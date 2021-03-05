//! VCF header contig record and key.

mod key;

use std::{convert::TryFrom, error, fmt, num};

use indexmap::IndexMap;

use super::{record, Record};

use self::key::Key;

/// A VCF header contig record (`contig`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Contig {
    id: String,
    len: Option<i32>,
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
    /// let contig = Contig::new(String::from("sq0"));
    /// ```
    pub fn new(id: String) -> Self {
        Self {
            id,
            len: None,
            fields: IndexMap::new(),
        }
    }

    /// Returns the ID of the contig.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Contig;
    /// let contig = Contig::new(String::from("sq0"));
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
    /// let contig = Contig::new(String::from("sq0"));
    /// assert_eq!(contig.len(), None);
    /// ```
    pub fn len(&self) -> Option<i32> {
        self.len
    }

    /// Returns the value of the field with the given key.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::convert::TryFrom;
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

        write!(f, "{}={}", Key::Id, self.id)?;

        if let Some(len) = self.len {
            write!(f, ",{}={}", Key::Length, len)?;
        }

        for (key, value) in &self.fields {
            write!(f, r#",{}="{}""#, key, value)?;
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
    /// A key is invalid.
    InvalidKey(key::ParseError),
    /// The length is invalid.
    InvalidLength(num::ParseIntError),
    /// A required field is missing.
    MissingField(Key),
}

impl error::Error for TryFromRecordError {}

impl fmt::Display for TryFromRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord => f.write_str("invalid record"),
            Self::MissingField(key) => write!(f, "missing {} field", key),
            Self::InvalidKey(e) => write!(f, "invalid key: {}", e),
            Self::InvalidLength(e) => write!(f, "invalid length: {}", e),
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
    let mut other_fields = IndexMap::new();

    for (raw_key, value) in fields {
        let key = raw_key.parse().map_err(TryFromRecordError::InvalidKey)?;

        match key {
            Key::Id => {
                id = Some(value);
            }
            Key::Length => {
                len = value
                    .parse()
                    .map(Some)
                    .map_err(TryFromRecordError::InvalidLength)?;
            }
            Key::Other(k) => {
                other_fields.insert(k, value);
            }
        }
    }

    Ok(Contig {
        id: id.ok_or(TryFromRecordError::MissingField(Key::Id))?,
        len,
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
                fields: vec![(
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
    fn test_try_from_record_for_contig_with_an_invalid_key() {
        let record = Record::new(
            record::Key::Contig,
            record::Value::Struct(vec![(String::new(), String::from("sq0"))]),
        );

        assert!(matches!(
            Contig::try_from(record),
            Err(TryFromRecordError::InvalidKey(_))
        ));
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
}
