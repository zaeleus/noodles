//! SAM header reference sequence builder.

use std::num::NonZero;

use super::ReferenceSequence;
use crate::header::record::value::map::{self, builder::BuildError};

/// A SAM header reference sequence builder.
#[derive(Debug, Default)]
pub struct Builder {
    length: Option<NonZero<usize>>,
}

impl map::Builder<ReferenceSequence> {
    /// Sets a reference sequence length.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZero;
    ///
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    ///
    /// let length = NonZero::try_from(13)?;
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
    ///     .set_length(length)
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.length(), length);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_length(mut self, length: NonZero<usize>) -> Self {
        self.inner.length = Some(length);
        self
    }
}

impl map::builder::Inner<ReferenceSequence> for Builder {
    fn build(self) -> Result<ReferenceSequence, BuildError> {
        let length = self.length.ok_or(BuildError::MissingField("LN"))?;

        Ok(ReferenceSequence { length })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();
        assert!(builder.length.is_none());
    }
}
