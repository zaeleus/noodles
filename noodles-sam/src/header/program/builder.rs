//! SAM header program and fields.

use std::{collections::HashMap, error, fmt};

use super::{Program, Tag};

/// A SAM header program builder.
#[derive(Debug, Default)]
pub struct Builder {
    id: Option<String>,
    name: Option<String>,
    command_line: Option<String>,
    previous_id: Option<String>,
    description: Option<String>,
    version: Option<String>,
    fields: HashMap<Tag, String>,
}

/// An error returned when a SAM header program fails to build.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BuildError {
    /// The ID is missing.
    MissingId,
}

impl error::Error for BuildError {}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingId => f.write_str("missing ID"),
        }
    }
}

impl Builder {
    /// Sets a program ID.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::program::builder;
    /// use noodles_sam::header::Program;
    /// let program = Program::builder().set_id("pg0").build()?;
    /// assert_eq!(program.id(), "pg0");
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_id<I>(mut self, id: I) -> Self
    where
        I: Into<String>,
    {
        self.id = Some(id.into());
        self
    }

    /// Sets a program name.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::program::builder;
    /// use noodles_sam::header::Program;
    ///
    /// let program = Program::builder()
    ///     .set_id("pg0")
    ///     .set_name("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(program.name(), Some("noodles"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_name<I>(mut self, name: I) -> Self
    where
        I: Into<String>,
    {
        self.name = Some(name.into());
        self
    }

    /// Sets a command line.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::program::builder;
    /// use noodles_sam::header::Program;
    ///
    /// let program = Program::builder()
    ///     .set_id("pg0")
    ///     .set_command_line("cargo run")
    ///     .build()?;
    ///
    /// assert_eq!(program.command_line(), Some("cargo run"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_command_line<I>(mut self, command_line: I) -> Self
    where
        I: Into<String>,
    {
        self.command_line = Some(command_line.into());
        self
    }

    /// Sets a previous program ID.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::program::builder;
    /// use noodles_sam::header::Program;
    ///
    /// let program = Program::builder()
    ///     .set_id("pg1")
    ///     .set_previous_id("pg0")
    ///     .build()?;
    ///
    /// assert_eq!(program.previous_id(), Some("pg0"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_previous_id<I>(mut self, previous_id: I) -> Self
    where
        I: Into<String>,
    {
        self.previous_id = Some(previous_id.into());
        self
    }

    /// Sets a description.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::program::builder;
    /// use noodles_sam::header::Program;
    ///
    /// let program = Program::builder()
    ///     .set_id("pg0")
    ///     .set_description("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(program.description(), Some("noodles"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_description<I>(mut self, description: I) -> Self
    where
        I: Into<String>,
    {
        self.description = Some(description.into());
        self
    }

    /// Sets a program version.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::program::builder;
    /// use noodles_sam::header::Program;
    ///
    /// let program = Program::builder()
    ///     .set_id("pg0")
    ///     .set_version("0.1.0")
    ///     .build()?;
    ///
    /// assert_eq!(program.version(), Some("0.1.0"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_version<I>(mut self, version: I) -> Self
    where
        I: Into<String>,
    {
        self.version = Some(version.into());
        self
    }

    /// Inserts a tag-raw value pair.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::program::builder;
    /// use noodles_sam::header::{program::Tag, Program};
    ///
    /// let zn = Tag::Other([b'z', b'n']);
    ///
    /// let program = Program::builder()
    ///     .set_id("pg0")
    ///     .insert(zn, "noodles")
    ///     .build()?;
    ///
    /// assert_eq!(program.fields().get(&zn), Some(&String::from("noodles")));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn insert<I>(mut self, tag: Tag, value: I) -> Self
    where
        I: Into<String>,
    {
        self.fields.insert(tag, value.into());
        self
    }

    /// Builds a SAM header program.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::program::builder;
    /// use noodles_sam::header::Program;
    /// let program = Program::builder().set_id("rg0").build()?;
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn build(self) -> Result<Program, BuildError> {
        let id = self.id.ok_or(BuildError::MissingId)?;

        Ok(Program {
            id,
            name: self.name,
            command_line: self.command_line,
            previous_id: self.previous_id,
            description: self.description,
            version: self.version,
            fields: self.fields,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert!(builder.id.is_none());
        assert!(builder.name.is_none());
        assert!(builder.command_line.is_none());
        assert!(builder.previous_id.is_none());
        assert!(builder.description.is_none());
        assert!(builder.version.is_none());
        assert!(builder.fields.is_empty());
    }

    #[test]
    fn test_build() {
        assert_eq!(Builder::default().build(), Err(BuildError::MissingId));
    }
}
