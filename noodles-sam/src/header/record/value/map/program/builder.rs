use super::Program;
use crate::header::record::value::map::{self, builder::BuildError};

/// A SAM header program builder.
#[derive(Debug, Default)]
pub struct Builder {
    name: Option<Vec<u8>>,
    command_line: Option<Vec<u8>>,
    previous_id: Option<Vec<u8>>,
    description: Option<Vec<u8>>,
    version: Option<Vec<u8>>,
}

impl map::Builder<Program> {
    /// Sets a program name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::record::value::{map::Program, Map}};
    ///
    /// let program = Map::<Program>::builder()
    ///     .set_name("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(program.name(), Some(&b"noodles"[..]));
    /// Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_name<I>(mut self, name: I) -> Self
    where
        I: Into<Vec<u8>>,
    {
        self.inner.name = Some(name.into());
        self
    }

    /// Sets a command line.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::record::value::{map::Program, Map}};
    ///
    /// let program = Map::<Program>::builder()
    ///     .set_command_line("cargo run")
    ///     .build()?;
    ///
    /// assert_eq!(program.command_line(), Some(&b"cargo run"[..]));
    /// Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_command_line<I>(mut self, command_line: I) -> Self
    where
        I: Into<Vec<u8>>,
    {
        self.inner.command_line = Some(command_line.into());
        self
    }

    /// Sets a previous program ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::record::value::{map::Program, Map}};
    ///
    /// let program = Map::<Program>::builder()
    ///     .set_previous_id("pg0")
    ///     .build()?;
    ///
    /// assert_eq!(program.previous_id(), Some(&b"pg0"[..]));
    /// Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_previous_id<I>(mut self, previous_id: I) -> Self
    where
        I: Into<Vec<u8>>,
    {
        self.inner.previous_id = Some(previous_id.into());
        self
    }

    /// Sets a description.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::record::value::{map::Program, Map}};
    ///
    /// let program = Map::<Program>::builder()
    ///     .set_description("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(program.description(), Some(&b"noodles"[..]));
    /// Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_description<I>(mut self, description: I) -> Self
    where
        I: Into<Vec<u8>>,
    {
        self.inner.description = Some(description.into());
        self
    }

    /// Sets a program version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::record::value::{map::Program, Map}};
    ///
    /// let program = Map::<Program>::builder()
    ///     .set_version("0.1.0")
    ///     .build()?;
    ///
    /// assert_eq!(program.version(), Some(&b"0.1.0"[..]));
    /// Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_version<I>(mut self, version: I) -> Self
    where
        I: Into<Vec<u8>>,
    {
        self.inner.version = Some(version.into());
        self
    }
}

impl map::builder::Inner<Program> for Builder {
    fn build(self) -> Result<Program, BuildError> {
        Ok(Program {
            name: self.name,
            command_line: self.command_line,
            previous_id: self.previous_id,
            description: self.description,
            version: self.version,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert!(builder.name.is_none());
        assert!(builder.command_line.is_none());
        assert!(builder.previous_id.is_none());
        assert!(builder.description.is_none());
        assert!(builder.version.is_none());
    }
}
