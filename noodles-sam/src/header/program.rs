//! SAM header program and fields.

mod tag;

use std::{collections::HashMap, convert::TryFrom, error, fmt};

pub use self::tag::Tag;

use super::record;

/// A SAM header program.
///
/// A program describes any program that created, viewed, or mutated a SAM file. The program ID is
/// guaranteed to be set.
#[derive(Debug)]
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
            ..Default::default()
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

    pub fn id_mut(&mut self) -> &mut String {
        &mut self.id
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

impl Default for Program {
    fn default() -> Self {
        Self {
            id: String::new(),
            fields: HashMap::new(),
        }
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
#[derive(Debug)]
pub enum ParseError {
    MissingRequiredTag(Tag),
    InvalidTag(tag::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingRequiredTag(tag) => write!(f, "missing required tag: {:?}", tag),
            Self::InvalidTag(e) => write!(f, "{}", e),
        }
    }
}

impl TryFrom<&[(String, String)]> for Program {
    type Error = ParseError;

    fn try_from(raw_fields: &[(String, String)]) -> Result<Self, Self::Error> {
        let mut program = Program::default();

        let mut has_id = false;

        for (raw_tag, value) in raw_fields {
            let tag = raw_tag.parse().map_err(ParseError::InvalidTag)?;

            if let Tag::Id = tag {
                program.id = value.into();
                has_id = true;
                continue;
            }

            program.fields.insert(tag, value.into());
        }

        if !has_id {
            return Err(ParseError::MissingRequiredTag(Tag::Id));
        }

        Ok(program)
    }
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
    fn test_from_str_with_no_id() {
        let fields = [(String::from("DS"), String::from("noodles"))];
        assert!(Program::try_from(&fields[..]).is_err());
    }
}
