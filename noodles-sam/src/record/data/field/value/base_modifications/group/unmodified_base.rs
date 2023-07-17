use std::{error, fmt};

/// The unmodified base as reported by the sequencer.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum UnmodifiedBase {
    /// Adenine.
    A,
    /// Cytosine.
    C,
    /// Guanine.
    G,
    /// Thymine.
    T,
    /// Uracil.
    U,
    /// Any base.
    N,
}

impl UnmodifiedBase {
    /// Complements the base.
    pub fn complement(&self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
            Self::U => Self::A,
            Self::N => Self::N,
        }
    }
}

/// An error returned when a base modifications group unmodified base fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => write!(f, "invalid input"),
        }
    }
}

impl TryFrom<u8> for UnmodifiedBase {
    type Error = ParseError;

    fn try_from(b: u8) -> Result<Self, Self::Error> {
        match b {
            b'A' => Ok(Self::A),
            b'C' => Ok(Self::C),
            b'G' => Ok(Self::G),
            b'T' => Ok(Self::T),
            b'U' => Ok(Self::U),
            b'N' => Ok(Self::N),
            _ => Err(ParseError::Invalid),
        }
    }
}

impl From<UnmodifiedBase> for crate::record::sequence::Base {
    fn from(unmodified_base: UnmodifiedBase) -> Self {
        match unmodified_base {
            UnmodifiedBase::A => Self::A,
            UnmodifiedBase::C => Self::C,
            UnmodifiedBase::G => Self::G,
            UnmodifiedBase::T => Self::T,
            UnmodifiedBase::U => Self::U,
            UnmodifiedBase::N => Self::N,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement() {
        assert_eq!(UnmodifiedBase::A.complement(), UnmodifiedBase::T);
        assert_eq!(UnmodifiedBase::C.complement(), UnmodifiedBase::G);
        assert_eq!(UnmodifiedBase::G.complement(), UnmodifiedBase::C);
        assert_eq!(UnmodifiedBase::T.complement(), UnmodifiedBase::A);
        assert_eq!(UnmodifiedBase::U.complement(), UnmodifiedBase::A);
        assert_eq!(UnmodifiedBase::N.complement(), UnmodifiedBase::N);
    }

    #[test]
    fn test_try_from_u8_for_unmodified_base() {
        fn t(b: u8, expected: UnmodifiedBase) {
            assert_eq!(UnmodifiedBase::try_from(b), Ok(expected));
        }

        t(b'A', UnmodifiedBase::A);
        t(b'C', UnmodifiedBase::C);
        t(b'G', UnmodifiedBase::G);
        t(b'T', UnmodifiedBase::T);
        t(b'U', UnmodifiedBase::U);
        t(b'N', UnmodifiedBase::N);

        assert_eq!(UnmodifiedBase::try_from(b'n'), Err(ParseError::Invalid));
    }

    #[test]
    fn test_from_unmodified_base_for_crate_sam_record_sequence_base() {
        use crate::record::sequence::Base;

        assert_eq!(Base::from(UnmodifiedBase::A), Base::A);
        assert_eq!(Base::from(UnmodifiedBase::C), Base::C);
        assert_eq!(Base::from(UnmodifiedBase::G), Base::G);
        assert_eq!(Base::from(UnmodifiedBase::T), Base::T);
        assert_eq!(Base::from(UnmodifiedBase::U), Base::U);
        assert_eq!(Base::from(UnmodifiedBase::N), Base::N);
    }
}
