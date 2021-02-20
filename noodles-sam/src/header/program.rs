//! SAM header program and fields.

mod builder;
mod tag;

use std::{collections::HashMap, convert::TryFrom, error, fmt};

pub use self::{builder::Builder, tag::Tag};

use super::{record, Record};

/// A SAM header program.
///
/// A program describes any program that created, viewed, or mutated a SAM file. The program ID is
/// guaranteed to be set.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Program {
    id: String,
    name: Option<String>,
    command_line: Option<String>,
    previous_id: Option<String>,
    description: Option<String>,
    version: Option<String>,
    fields: HashMap<Tag, String>,
}

impl Program {
    /// Creates a SAM header program builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::Program;
    /// let builder = Program::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

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
            name: None,
            command_line: None,
            previous_id: None,
            description: None,
            version: None,
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

    /// Returns the program name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::Program;
    /// let program = Program::new(String::from("pg0"));
    /// assert!(program.name().is_none());
    /// ```
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    /// Returns the command line.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::Program;
    /// let program = Program::new(String::from("pg0"));
    /// assert!(program.command_line().is_none());
    /// ```
    pub fn command_line(&self) -> Option<&str> {
        self.command_line.as_deref()
    }

    /// Returns the previous program ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::Program;
    /// let program = Program::new(String::from("pg0"));
    /// assert!(program.previous_id().is_none());
    /// ```
    pub fn previous_id(&self) -> Option<&str> {
        self.previous_id.as_deref()
    }

    /// Returns the description.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::Program;
    /// let program = Program::new(String::from("pg0"));
    /// assert!(program.description().is_none());
    /// ```
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }

    /// Returns the program version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::Program;
    /// let program = Program::new(String::from("pg0"));
    /// assert!(program.version().is_none());
    /// ```
    pub fn version(&self) -> Option<&str> {
        self.version.as_deref()
    }

    /// Returns the raw fields of the program.
    ///
    /// This includes any field that is not specially handled by the structure itself. For example,
    /// this will not include the ID field, as it is parsed and available as [`Self::id`].
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{program::Tag, Program};
    ///
    /// let program = Program::builder()
    ///     .set_id(String::from("pg0"))
    ///     .insert(Tag::Other(String::from("zn")), String::from("noodles"))
    ///     .build();
    ///
    /// let fields = program.fields();
    /// assert_eq!(fields.len(), 1);
    /// assert_eq!(
    ///     fields.get(&Tag::Other(String::from("zn"))),
    ///     Some(&String::from("noodles"))
    /// );
    ///
    /// assert_eq!(fields.get(&Tag::Id), None);
    /// assert_eq!(program.id(), "pg0");
    /// ```
    pub fn fields(&self) -> &HashMap<Tag, String> {
        &self.fields
    }
}

impl fmt::Display for Program {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", record::Kind::Program)?;
        write!(f, "\t{}:{}", Tag::Id, self.id)?;

        if let Some(name) = self.name() {
            write!(f, "\t{}:{}", Tag::Name, name)?;
        }

        if let Some(command_line) = self.command_line() {
            write!(f, "\t{}:{}", Tag::CommandLine, command_line)?;
        }

        if let Some(previous_id) = self.previous_id() {
            write!(f, "\t{}:{}", Tag::PreviousId, previous_id)?;
        }

        if let Some(description) = self.description() {
            write!(f, "\t{}:{}", Tag::Description, description)?;
        }

        if let Some(version) = self.version() {
            write!(f, "\t{}:{}", Tag::Version, version)?;
        }

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
            Self::InvalidTag(e) => write!(f, "invalid tag: {}", e),
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
    let mut builder = Program::builder();
    let mut id = None;

    for (raw_tag, value) in raw_fields {
        let tag = raw_tag.parse().map_err(TryFromRecordError::InvalidTag)?;

        builder = match tag {
            Tag::Id => {
                id = Some(value);
                builder
            }
            Tag::Name => builder.set_name(value),
            Tag::CommandLine => builder.set_command_line(value),
            Tag::PreviousId => builder.set_previous_id(value),
            Tag::Description => builder.set_description(value),
            Tag::Version => builder.set_version(value),
            _ => builder.insert(tag, value),
        }
    }

    if let Some(i) = id {
        builder = builder.set_id(i);
    } else {
        return Err(TryFromRecordError::MissingRequiredTag(Tag::Id));
    }

    Ok(builder.build())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let program = Program::builder().set_id("pg0").set_name("noodles").build();
        assert_eq!(program.to_string(), "@PG\tID:pg0\tPN:noodles");
    }

    #[test]
    fn test_try_from_record_for_program_with_invalid_record() {
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
    fn test_try_from_record_for_program_with_no_id() {
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
