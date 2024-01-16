use super::{Header, Version};
use crate::header::record::value::map::{self, builder::BuildError};

/// A SAM header header builder.
#[derive(Debug, Default)]
pub struct Builder {
    version: Option<Version>,
}

impl map::Builder<Header> {
    /// Sets a format version.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{
    ///     map::{self, header::Version},
    ///     Map,
    /// };
    ///
    /// let version = Version::new(1, 6);
    /// let header = Map::<map::Header>::builder().set_version(version).build()?;
    /// assert_eq!(header.version(), version);
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_version(mut self, version: Version) -> Self {
        self.inner.version = Some(version);
        self
    }
}

impl map::builder::Inner<Header> for Builder {
    fn build(self) -> Result<Header, BuildError> {
        Ok(Header {
            version: self.version.unwrap_or_default(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();
        assert!(builder.version.is_none());
    }
}
