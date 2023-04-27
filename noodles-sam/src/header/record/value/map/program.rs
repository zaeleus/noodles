//! SAM header record program map value.

mod builder;
mod tag;

use std::fmt;

use self::{
    builder::Builder,
    tag::{StandardTag, Tag},
};
use super::{Fields, Inner, Map, TryFromFieldsError};

// A SAM header record program map value.
///
/// A program describes any program that created, viewed, or mutated a SAM file.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Program {
    name: Option<String>,
    command_line: Option<String>,
    previous_id: Option<String>,
    description: Option<String>,
    version: Option<String>,
}

impl Inner for Program {
    type StandardTag = StandardTag;
    type Builder = Builder;
}

impl Map<Program> {
    /// Returns the program name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::Program, Map};
    /// let program = Map::<Program>::default();
    /// assert!(program.name().is_none());
    /// ```
    pub fn name(&self) -> Option<&str> {
        self.inner.name.as_deref()
    }

    /// Returns the command line.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::Program, Map};
    /// let program = Map::<Program>::default();
    /// assert!(program.command_line().is_none());
    /// ```
    pub fn command_line(&self) -> Option<&str> {
        self.inner.command_line.as_deref()
    }

    /// Returns the previous program ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::Program, Map};
    /// let program = Map::<Program>::default();
    /// assert!(program.previous_id().is_none());
    /// ```
    pub fn previous_id(&self) -> Option<&str> {
        self.inner.previous_id.as_deref()
    }

    /// Returns the description.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::Program, Map};
    /// let program = Map::<Program>::default();
    /// assert!(program.description().is_none());
    /// ```
    pub fn description(&self) -> Option<&str> {
        self.inner.description.as_deref()
    }

    /// Returns the program version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::Program, Map};
    /// let program = Map::<Program>::default();
    /// assert!(program.version().is_none());
    /// ```
    pub fn version(&self) -> Option<&str> {
        self.inner.version.as_deref()
    }
}

impl fmt::Display for Map<Program> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(name) = self.name() {
            write!(f, "\t{}:{name}", tag::NAME)?;
        }

        if let Some(command_line) = self.command_line() {
            write!(f, "\t{}:{command_line}", tag::COMMAND_LINE)?;
        }

        if let Some(previous_id) = self.previous_id() {
            write!(f, "\t{}:{previous_id}", tag::PREVIOUS_ID)?;
        }

        if let Some(description) = self.description() {
            write!(f, "\t{}:{description}", tag::DESCRIPTION)?;
        }

        if let Some(version) = self.version() {
            write!(f, "\t{}:{version}", tag::VERSION)?;
        }

        super::fmt_display_other_fields(f, self.other_fields())?;

        Ok(())
    }
}

impl TryFrom<Fields> for Map<Program> {
    type Error = TryFromFieldsError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        let mut name = None;
        let mut command_line = None;
        let mut previous_id = None;
        let mut description = None;
        let mut version = None;

        let mut other_fields = super::init_other_fields();

        for (key, value) in fields {
            let tag = key.parse().map_err(|_| TryFromFieldsError::InvalidTag)?;

            match tag {
                tag::ID => return Err(TryFromFieldsError::DuplicateTag),
                tag::NAME => name = Some(value),
                tag::COMMAND_LINE => command_line = Some(value),
                tag::PREVIOUS_ID => previous_id = Some(value),
                tag::DESCRIPTION => description = Some(value),
                tag::VERSION => version = Some(value),
                Tag::Other(t) => super::insert_other_field(&mut other_fields, t, value)?,
            }
        }

        Ok(Self {
            inner: Program {
                name,
                command_line,
                previous_id,
                description,
                version,
            },
            other_fields,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::header::record::value::map::builder::BuildError;

    #[test]
    fn test_fmt() -> Result<(), BuildError> {
        let program = Map::<Program>::builder()
            .set_name("noodles-sam")
            .set_version("0.23.0")
            .build()?;

        assert_eq!(program.to_string(), "\tPN:noodles-sam\tVN:0.23.0");

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_program_with_duplicate_id() {
        let fields = vec![(String::from("ID"), String::from("pg0"))];

        assert_eq!(
            Map::<Program>::try_from(fields),
            Err(TryFromFieldsError::DuplicateTag)
        );
    }
}
