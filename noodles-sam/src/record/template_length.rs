//! SAM record template length.

use std::{
    error, fmt,
    num::{self, NonZeroUsize},
    ops::Not,
};

/// A SAM record template length.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum TemplateLength {
    /// The leftmost segment.
    Left(NonZeroUsize),
    /// The rightmost segment.
    Right(NonZeroUsize),
}

impl TemplateLength {
    /// Creates a template length for the leftmost segment.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::TemplateLength;
    ///
    /// let template_length = TemplateLength::left(0);
    /// assert!(template_length.is_none());
    ///
    /// let template_length = TemplateLength::left(8);
    /// assert!(matches!(template_length, Some(TemplateLength::Left(_))));
    /// assert_eq!(template_length.map(usize::from), Some(8));
    /// ```
    pub const fn left(len: usize) -> Option<Self> {
        if let Some(n) = NonZeroUsize::new(len) {
            Some(Self::Left(n))
        } else {
            None
        }
    }

    /// Creates a template length for the rightmost segment.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::TemplateLength;
    ///
    /// let template_length = TemplateLength::right(0);
    /// assert!(template_length.is_none());
    ///
    /// let template_length = TemplateLength::right(8);
    /// assert!(matches!(template_length, Some(TemplateLength::Right(_))));
    /// assert_eq!(template_length.map(usize::from), Some(8));
    /// ```
    pub const fn right(len: usize) -> Option<Self> {
        if let Some(n) = NonZeroUsize::new(len) {
            Some(Self::Right(n))
        } else {
            None
        }
    }
}

impl fmt::Display for TemplateLength {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Left(len) => len.fmt(f),
            Self::Right(len) => write!(f, "-{}", len),
        }
    }
}

impl Not for TemplateLength {
    type Output = Self;

    fn not(self) -> Self::Output {
        match self {
            Self::Left(len) => Self::Right(len),
            Self::Right(len) => Self::Left(len),
        }
    }
}

/// An error returned when a raw SAM record template length fails to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromIntError {
    /// The input is invalid.
    Invalid(num::TryFromIntError),
}

impl error::Error for TryFromIntError {}

impl fmt::Display for TryFromIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(e) => write!(f, "invalid input: {}", e),
        }
    }
}

impl TryFrom<i32> for TemplateLength {
    type Error = TryFromIntError;

    fn try_from(n: i32) -> Result<Self, Self::Error> {
        usize::try_from(abs(n))
            .and_then(NonZeroUsize::try_from)
            .map(|len| {
                if n > 0 {
                    Self::Left(len)
                } else {
                    Self::Right(len)
                }
            })
            .map_err(TryFromIntError::Invalid)
    }
}

impl From<TemplateLength> for usize {
    fn from(template_length: TemplateLength) -> Self {
        match template_length {
            TemplateLength::Left(len) => usize::from(len),
            TemplateLength::Right(len) => usize::from(len),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_not() -> Result<(), TryFromIntError> {
        let left_template_length = TemplateLength::try_from(8)?;
        let right_template_length = TemplateLength::try_from(-8)?;

        assert_eq!(!left_template_length, right_template_length);
        assert_eq!(!right_template_length, left_template_length);

        Ok(())
    }

    #[test]
    fn test_fmt() -> Result<(), TryFromIntError> {
        let template_length = TemplateLength::try_from(8)?;
        assert_eq!(template_length.to_string(), "8");

        let template_length = TemplateLength::try_from(-8)?;
        assert_eq!(template_length.to_string(), "-8");

        Ok(())
    }

    #[test]
    fn test_try_from_i32_for_template_length() -> Result<(), num::TryFromIntError> {
        assert_eq!(
            TemplateLength::try_from(8),
            Ok(TemplateLength::Left(NonZeroUsize::try_from(8)?))
        );

        assert_eq!(
            TemplateLength::try_from(-8),
            Ok(TemplateLength::Right(NonZeroUsize::try_from(8)?))
        );

        assert_eq!(
            TemplateLength::try_from(i32::MAX),
            Ok(TemplateLength::Left(NonZeroUsize::try_from((1 << 31) - 1)?))
        );

        assert_eq!(
            TemplateLength::try_from(i32::MIN),
            Ok(TemplateLength::Right(NonZeroUsize::try_from(1 << 31)?))
        );

        assert!(matches!(
            TemplateLength::try_from(0),
            Err(TryFromIntError::Invalid(_))
        ));

        Ok(())
    }

    #[test]
    fn test_from_template_length_for_usize() -> Result<(), TryFromIntError> {
        assert_eq!(usize::from(TemplateLength::try_from(8)?), 8);
        assert_eq!(usize::from(TemplateLength::try_from(-8)?), 8);
        Ok(())
    }
}

// Returns the absolute value without overflow.
fn abs(n: i32) -> u32 {
    if n < 0 {
        0u32.wrapping_sub(n as u32)
    } else {
        n as u32
    }
}
