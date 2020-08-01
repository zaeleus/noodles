//! SAM header program and fields.

mod tag;

use std::{collections::HashMap, convert::TryFrom, error, fmt};

pub use self::tag::Tag;

use super::{record, Record};

/// A SAM header program.
///
/// A program describes any program that created, viewed, or mutated a SAM file. The program ID is
/// guaranteed to be set.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Program {
    id: String,
    fields: HashMap<Tag, String>,
}

impl Program {
    /// Creates a program with an ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::Program;
    /// let program = Program::new(String::from("pg0"));
    /// assert_eq!(program.id(), "pg0");
    /// ```
    pub fn new(id: String) -> Self {
        Self {
            id,
            fields: HashMap::new(),
        }
    }

    /// Returns the program ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::Program;
    /// let program = Program::new(String::from("pg0"));
    /// assert_eq!(program.id(), "pg0");
    /// ```
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Returns a mutable reference to the program ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::Program;
    ///
    /// let mut program = Program::new(String::from("pg0"));
    /// assert_eq!(program.id(), "pg0");
    ///
    /// *program.id_mut() = String::from("pg1");
    /// assert_eq!(program.id(), "pg1");
    /// ```
    pub fn id_mut(&mut self) -> &mut String {
        &mut self.id
    }

    /// Returns the raw fields of the program.
    ///
    /// This includes any field that is not specially handled by the structure itself. For example,
    /// this will not include the ID field, as it is parsed and available as [`id`].
    ///
    /// [`id`]: #method.id
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{program::Tag, Program};
    ///
    /// let mut program = Program::new(String::from("pg0"));
    /// program.insert(Tag::Name, String::from("noodles-sam"));
    ///
    /// let fields = program.fields();
    /// assert_eq!(fields.len(), 1);
    /// assert_eq!(fields.get(&Tag::Name), Some(&String::from("noodles-sam")));
    ///
    /// assert_eq!(fields.get(&Tag::Id), None);
    /// assert_eq!(program.id(), "pg0");
    /// ```
    pub fn fields(&self) -> &HashMap<Tag, String> {
        &self.fields
    }

    /// Returns a reference to the raw field value mapped to the given key.
    ///
    /// This can only be used for fields with unparsed values. For a program, [`id`] must be used
    /// instead of `get(Tag::Id)`.
    ///
    /// [`id`]: #method.id
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{program::Tag, Program};
    ///
    /// let mut program = Program::new(String::from("pg0"));
    /// program.insert(Tag::Name, String::from("noodles-sam"));
    ///
    /// assert_eq!(program.get(&Tag::Name), Some(&String::from("noodles-sam")));
    /// assert_eq!(program.get(&Tag::Id), None);
    /// ```
    pub fn get(&self, tag: &Tag) -> Option<&String> {
        self.fields.get(tag)
    }

    /// Inserts a tag-raw value pair into the program.
    ///
    /// This follows similar semantics to [`std::collections::HashMap::insert`].
    ///
    /// [`std::collections::HashMap::insert`]: https://doc.rust-lang.org/stable/std/collections/struct.HashMap.html#method.insert
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{program::Tag, Program};
    /// let mut program = Program::new(String::from("pg0"));
    /// program.insert(Tag::Name, String::from("noodles-sam"));
    /// ```
    pub fn insert(&mut self, tag: Tag, value: String) -> Option<String> {
        self.fields.insert(tag, value)
    }
}

impl fmt::Display for Program {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", record::Kind::Program)?;
        write!(f, "\t{}:{}", Tag::Id, self.id)?;

        for (tag, value) in &self.fields {
            write!(f, "\t{}:{}", tag, value)?;
        }

        Ok(())
    }
}

/// An error returned when a raw SAM header program fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromRecordError {
    /// The record is invalid.
    InvalidRecord,
    /// A required tag is missing.
    MissingRequiredTag(Tag),
    /// A tag is invalid.
    InvalidTag(tag::ParseError),
}

impl error::Error for TryFromRecordError {}

impl fmt::Display for TryFromRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord => f.write_str("invalid record"),
            Self::MissingRequiredTag(tag) => write!(f, "missing required tag: {:?}", tag),
            Self::InvalidTag(e) => write!(f, "{}", e),
        }
    }
}

impl TryFrom<Record> for Program {
    type Error = TryFromRecordError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        match record.into() {
            (record::Kind::Program, record::Value::Map(fields)) => parse_map(fields),
            _ => Err(TryFromRecordError::InvalidRecord),
        }
    }
}

fn parse_map(raw_fields: Vec<(String, String)>) -> Result<Program, TryFromRecordError> {
    let mut id = None;
    let mut fields = HashMap::new();

    for (raw_tag, value) in raw_fields {
        let tag = raw_tag.parse().map_err(TryFromRecordError::InvalidTag)?;

        if let Tag::Id = tag {
            id = Some(value);
        } else {
            fields.insert(tag, value);
        }
    }

    Ok(Program {
        id: id.ok_or_else(|| TryFromRecordError::MissingRequiredTag(Tag::Id))?,
        fields,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let mut program = Program::new(String::from("pg0"));

        program.fields.insert(Tag::Name, String::from("noodles"));

        let actual = format!("{}", program);
        let expected = "@PG\tID:pg0\tPN:noodles";

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_from_str_with_invalid_record() {
        let record = Record::new(
            record::Kind::Comment,
            record::Value::String(String::from("noodles-sam")),
        );

        assert_eq!(
            Program::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_from_str_with_no_id() {
        let record = Record::new(
            record::Kind::Program,
            record::Value::Map(vec![(String::from("PN"), String::from("noodles"))]),
        );

        assert_eq!(
            Program::try_from(record),
            Err(TryFromRecordError::MissingRequiredTag(Tag::Id))
        );
    }
}
