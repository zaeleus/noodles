use std::{convert::TryFrom, error, fmt};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Base {
    A,
    C,
    G,
    T,
    N,
}

#[derive(Debug)]
pub struct TryFromCharError(char);

impl error::Error for TryFromCharError {}

impl fmt::Display for TryFromCharError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "invalid reference base: expected {{A, C, G, T, N}}, got {}",
            self.0
        )
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
    fn test_try_from_char_for_base() -> Result<(), TryFromCharError> {
        assert_eq!(Base::try_from('A')?, Base::A);
        assert_eq!(Base::try_from('C')?, Base::C);
        assert_eq!(Base::try_from('G')?, Base::G);
        assert_eq!(Base::try_from('T')?, Base::T);
        assert_eq!(Base::try_from('N')?, Base::N);

        assert!(Base::try_from('Z').is_err());

        Ok(())
    }
}
