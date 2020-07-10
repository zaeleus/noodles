use std::{convert::TryFrom, error, fmt};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Base {
    A,
    C,
    G,
    T,
    N,
}

impl From<Base> for char {
    fn from(base: Base) -> Self {
        match base {
            Base::A => 'A',
            Base::C => 'C',
            Base::G => 'G',
            Base::T => 'T',
            Base::N => 'N',
        }
    }
}

/// An error returned when a character fails to convert to a base.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromCharError(char);

impl error::Error for TryFromCharError {}

impl fmt::Display for TryFromCharError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "expected {{A, C, G, T, N}}, got {}", self.0)
    }
}

impl TryFrom<char> for Base {
    type Error = TryFromCharError;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            'A' => Ok(Self::A),
            'C' => Ok(Self::C),
            'G' => Ok(Self::G),
            'T' => Ok(Self::T),
            'N' => Ok(Self::N),
            _ => Err(TryFromCharError(value)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_base_for_char() {
        assert_eq!(char::from(Base::A), 'A');
        assert_eq!(char::from(Base::C), 'C');
        assert_eq!(char::from(Base::G), 'G');
        assert_eq!(char::from(Base::T), 'T');
        assert_eq!(char::from(Base::N), 'N');
    }

    #[test]
    fn test_try_from_char_for_base() {
        assert_eq!(Base::try_from('A'), Ok(Base::A));
        assert_eq!(Base::try_from('C'), Ok(Base::C));
        assert_eq!(Base::try_from('G'), Ok(Base::G));
        assert_eq!(Base::try_from('T'), Ok(Base::T));
        assert_eq!(Base::try_from('N'), Ok(Base::N));

        assert_eq!(Base::try_from('Z'), Err(TryFromCharError('Z')));
    }
}
