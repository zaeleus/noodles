mod other;

use std::{error, fmt, marker::PhantomData, str::FromStr};

pub use self::other::Other;

pub(crate) const LENGTH: usize = 2;

pub trait Standard: AsRef<[u8; LENGTH]> + TryFrom<[u8; LENGTH]> {}

/// A SAM header record map value field tag.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Tag<S> {
    /// A standard tag.
    Standard(S),
    /// A nonstandard tag.
    Other(Other<S>),
}

impl<S> fmt::Display for Tag<S>
where
    S: Standard,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Standard(tag) => {
                let b = tag.as_ref();
                char::from(b[0]).fmt(f)?;
                char::from(b[1]).fmt(f)?;
            }
            Self::Other(tag) => {
                tag.fmt(f)?;
            }
        }

        Ok(())
    }
}

impl<S> From<[u8; LENGTH]> for Tag<S>
where
    S: Standard,
{
    fn from(s: [u8; LENGTH]) -> Self {
        match S::try_from(s) {
            Ok(tag) => Self::Standard(tag),
            Err(_) => Self::Other(Other(s, PhantomData)),
        }
    }
}

/// An error returned when a raw SAM header map tag fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is not two bytes.
    LengthMismatch { actual: usize },
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::LengthMismatch { actual } => {
                write!(f, "length mismatch: expected {LENGTH}, got {actual}")
            }
            Self::Invalid => write!(f, "invalid input"),
        }
    }
}

impl<S> FromStr for Tag<S>
where
    S: Standard,
{
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            Err(ParseError::Empty)
        } else if s.len() != 2 {
            Err(ParseError::LengthMismatch { actual: s.len() })
        } else {
            // SAFETY: `s` is 2 bytes.
            let raw_tag = s.as_bytes().try_into().unwrap();

            if is_valid_tag(raw_tag) {
                Ok(Self::from(raw_tag))
            } else {
                Err(ParseError::Invalid)
            }
        }
    }
}

fn is_valid_tag(b: [u8; LENGTH]) -> bool {
    b[0].is_ascii_alphabetic() && b[1].is_ascii_alphanumeric()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str_for_tag() {
        const ID: [u8; LENGTH] = [b'I', b'D'];

        #[derive(Debug, Eq, PartialEq)]
        enum TestTag {
            Id,
        }

        impl Standard for TestTag {}

        impl AsRef<[u8; LENGTH]> for TestTag {
            fn as_ref(&self) -> &[u8; LENGTH] {
                &ID
            }
        }

        impl TryFrom<[u8; LENGTH]> for TestTag {
            type Error = ();

            fn try_from(b: [u8; LENGTH]) -> Result<Self, Self::Error> {
                match b {
                    ID => Ok(Self::Id),
                    _ => Err(()),
                }
            }
        }

        assert_eq!("ID".parse(), Ok(Tag::Standard(TestTag::Id)));
        assert_eq!(
            "VN".parse(),
            Ok(Tag::<TestTag>::Other(Other([b'V', b'N'], PhantomData)))
        );

        assert_eq!(
            "V".parse::<Tag<TestTag>>(),
            Err(ParseError::LengthMismatch { actual: 1 })
        );
        assert_eq!("0V".parse::<Tag<TestTag>>(), Err(ParseError::Invalid));
        assert_eq!(
            "VER".parse::<Tag<TestTag>>(),
            Err(ParseError::LengthMismatch { actual: 3 })
        );
    }
}
