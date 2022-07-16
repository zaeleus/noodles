//! VCF header contig record and key.

mod builder;
pub mod name;
mod tag;

pub use self::{name::Name, tag::Tag};

use std::{error, fmt, num};

use indexmap::IndexMap;

use self::builder::Builder;
use super::record;

/// A VCF header contig record (`contig`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Contig {
    id: Name,
    len: Option<usize>,
    idx: Option<usize>,
    fields: IndexMap<tag::Other, String>,
}

#[allow(clippy::len_without_is_empty)]
impl Contig {
    pub(crate) fn try_from_fields(
        id: String,
        fields: IndexMap<String, String>,
    ) -> Result<Self, TryFromRecordError> {
        parse_struct(id, fields)
    }

    /// Creates a VCF header contig record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Contig;
    /// let contig = Contig::new("sq0".parse()?);
    /// Ok::<_, noodles_vcf::header::contig::name::ParseError>(())
    /// ```
    pub fn new(id: Name) -> Self {
        Self {
            id,
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
    /// use noodles_vcf::header::{contig::Name, Contig};
    /// let name: Name = "sq0".parse()?;
    /// let contig = Contig::new(name.clone());
    /// assert_eq!(contig.id(), &name);
    /// Ok::<_, noodles_vcf::header::contig::name::ParseError>(())
    /// ```
    pub fn id(&self) -> &Name {
        &self.id
    }

    /// Returns the length of the contig, if it is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Contig;
    /// let contig = Contig::new("sq0".parse()?);
    /// assert_eq!(contig.len(), None);
    /// Ok::<_, noodles_vcf::header::contig::name::ParseError>(())
    /// ```
    pub fn len(&self) -> Option<usize> {
        self.len
    }

    /// Returns a mutable reference to the length of the contig.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Contig;
    ///
    /// let mut contig = Contig::new("sq0".parse()?);
    /// assert!(contig.len().is_none());
    ///
    /// *contig.len_mut() = Some(8);
    /// assert_eq!(contig.len(), Some(8));
    /// Ok::<_, noodles_vcf::header::contig::name::ParseError>(())
    /// ```
    pub fn len_mut(&mut self) -> &mut Option<usize> {
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
    /// let contig = Contig::new("sq0".parse()?);
    /// assert!(contig.idx().is_none());
    /// Ok::<_, noodles_vcf::header::contig::name::ParseError>(())
    /// ```
    pub fn idx(&self) -> Option<usize> {
        self.idx
    }

    /// Returns the value of the field with the given key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::Record;
    ///
    /// let s = "##contig=<ID=sq0,md5=d7eba311421bbc9d3ada44709dd61534>";
    /// let record = s.parse()?;
    ///
    /// let contig = match record {
    ///     Record::Contig(contig) => contig,
    ///     _ => unreachable!(),
    /// };
    ///
    /// assert_eq!(contig.get("md5"), Some("d7eba311421bbc9d3ada44709dd61534"));
    /// assert!(contig.get("species").is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn get(&self, key: &str) -> Option<&str> {
        self.fields.get(key).map(|s| &**s)
    }
}

impl fmt::Display for Contig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(record::PREFIX)?;
        f.write_str(record::key::CONTIG.as_ref())?;
        f.write_str("=<")?;

        write!(f, "{}={}", tag::ID, self.id)?;

        if let Some(len) = self.len {
            write!(f, ",{}={}", tag::LENGTH, len)?;
        }

        for (key, value) in &self.fields {
            write!(f, ",{}=", key)?;
            super::fmt::write_escaped_string(f, value)?;
        }

        if let Some(idx) = self.idx() {
            write!(f, ",{}={}", tag::IDX, idx)?;
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
    InvalidId(name::ParseError),
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
            Self::InvalidId(e) => write!(f, "invalid ID: {}", e),
            Self::InvalidLength(e) => write!(f, "invalid length: {}", e),
            Self::InvalidIdx(e) => write!(f, "invalid index (`{}`): {}", tag::IDX, e),
        }
    }
}

fn parse_struct(
    raw_id: String,
    fields: IndexMap<String, String>,
) -> Result<Contig, TryFromRecordError> {
    let mut builder = Builder::default();

    let id = raw_id.parse().map_err(TryFromRecordError::InvalidId)?;
    builder = builder.set_id(id);

    for (key, value) in fields {
        let tag = Tag::from(key);

        builder = match tag {
            tag::ID => todo!(),
            tag::LENGTH => {
                let len = value.parse().map_err(TryFromRecordError::InvalidLength)?;
                builder.set_len(len)
            }
            tag::IDX => {
                let idx = value.parse().map_err(TryFromRecordError::InvalidIdx)?;
                builder.set_idx(idx)
            }
            Tag::Other(t) => builder.insert(t, value),
        }
    }

    builder
        .build()
        .map_err(|_| TryFromRecordError::InvalidRecord)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_record() -> (String, IndexMap<String, String>) {
        (
            String::from("sq0"),
            [
                (String::from("length"), String::from("13")),
                (
                    String::from("md5"),
                    String::from("d7eba311421bbc9d3ada44709dd61534"),
                ),
            ]
            .into_iter()
            .collect(),
        )
    }

    #[test]
    fn test_fmt() -> Result<(), TryFromRecordError> {
        let (id, fields) = build_record();
        let contig = Contig::try_from_fields(id, fields)?;

        let expected = r#"##contig=<ID=sq0,length=13,md5="d7eba311421bbc9d3ada44709dd61534">"#;
        assert_eq!(contig.to_string(), expected);

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_contig() -> Result<(), Box<dyn std::error::Error>> {
        let (id, fields) = build_record();

        assert_eq!(
            Contig::try_from_fields(id, fields),
            Ok(Contig {
                id: "sq0".parse()?,
                len: Some(13),
                idx: None,
                fields: [(
                    Tag::other("md5").ok_or("invalid tag")?,
                    String::from("d7eba311421bbc9d3ada44709dd61534")
                )]
                .into_iter()
                .collect(),
            })
        );

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_contig_with_extra_fields() -> Result<(), Box<dyn std::error::Error>>
    {
        let id = String::from("sq0");
        let fields = [
            (String::from("length"), String::from("13")),
            (
                String::from("md5"),
                String::from("d7eba311421bbc9d3ada44709dd61534"),
            ),
            (String::from("IDX"), String::from("1")),
        ]
        .into_iter()
        .collect();

        assert_eq!(
            Contig::try_from_fields(id, fields),
            Ok(Contig {
                id: "sq0".parse()?,
                len: Some(13),
                idx: Some(1),
                fields: [(
                    Tag::other("md5").ok_or("invalid tag")?,
                    String::from("d7eba311421bbc9d3ada44709dd61534")
                )]
                .into_iter()
                .collect(),
            })
        );

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_contig_with_an_invalid_id() {
        assert!(matches!(
            Contig::try_from_fields(String::from("sq 0"), IndexMap::new()),
            Err(TryFromRecordError::InvalidId(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_contig_with_an_invalid_length() {
        let id = String::from("sq0");
        let fields = [(String::from("length"), String::from("eight"))]
            .into_iter()
            .collect();

        assert!(matches!(
            Contig::try_from_fields(id, fields),
            Err(TryFromRecordError::InvalidLength(_))
        ));
    }

    #[test]
    fn test_try_from_record_for_contig_with_an_invalid_idx() {
        let id = String::from("sq0");
        let fields = [
            (String::from("length"), String::from("13")),
            (String::from("IDX"), String::from("zero")),
        ]
        .into_iter()
        .collect();

        assert!(matches!(
            Contig::try_from_fields(id, fields),
            Err(TryFromRecordError::InvalidIdx(_))
        ));
    }
}
