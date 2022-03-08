use std::{error, fmt};

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum Base {
    A,
    C,
    G,
    T,
    N,
}

impl Default for Base {
    fn default() -> Self {
        Self::N
    }
}

#[derive(Debug, Eq, PartialEq)]
pub struct TryFromIntError(u8);

impl fmt::Display for TryFromIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid substitution matrix base: expected {{A, C, G, T, N}}, got {:#04x}",
            self.0
        )
    }
}

impl error::Error for TryFromIntError {}

impl TryFrom<u8> for Base {
    type Error = TryFromIntError;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        match n {
            b'A' => Ok(Self::A),
            b'C' => Ok(Self::C),
            b'G' => Ok(Self::G),
            b'T' => Ok(Self::T),
            b'N' => Ok(Self::N),
            _ => Err(TryFromIntError(n)),
        }
    }
}

impl From<Base> for u8 {
    fn from(base: Base) -> Self {
        match base {
            Base::A => b'A',
            Base::C => b'C',
            Base::G => b'G',
            Base::T => b'T',
            Base::N => b'N',
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(Base::default(), Base::N);
    }

    #[test]
    fn test_try_from_u8_for_base() {
        assert_eq!(Base::try_from(b'A'), Ok(Base::A));
        assert_eq!(Base::try_from(b'C'), Ok(Base::C));
        assert_eq!(Base::try_from(b'G'), Ok(Base::G));
        assert_eq!(Base::try_from(b'T'), Ok(Base::T));
        assert_eq!(Base::try_from(b'N'), Ok(Base::N));
        assert_eq!(Base::try_from(b'U'), Err(TryFromIntError(b'U')));
    }

    #[test]
    fn test_from_base_for_char() {
        assert_eq!(u8::from(Base::A), b'A');
        assert_eq!(u8::from(Base::C), b'C');
        assert_eq!(u8::from(Base::G), b'G');
        assert_eq!(u8::from(Base::T), b'T');
        assert_eq!(u8::from(Base::N), b'N');
    }
}
