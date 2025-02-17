use std::io;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum Base {
    A,
    C,
    G,
    T,
    N,
}

impl TryFrom<u8> for Base {
    type Error = io::Error;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        match n.to_ascii_uppercase() {
            b'A' => Ok(Self::A),
            b'C' => Ok(Self::C),
            b'G' => Ok(Self::G),
            b'T' => Ok(Self::T),
            b'N' => Ok(Self::N),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid substitution base",
            )),
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
    fn test_try_from_u8_for_base() -> io::Result<()> {
        fn t(n: u8, expected: Base) -> io::Result<()> {
            assert_eq!(Base::try_from(n)?, expected);
            assert_eq!(Base::try_from(n.to_ascii_lowercase())?, expected);
            Ok(())
        }

        t(b'A', Base::A)?;
        t(b'C', Base::C)?;
        t(b'G', Base::G)?;
        t(b'T', Base::T)?;
        t(b'N', Base::N)?;

        assert!(matches!(
            Base::try_from(b'U'),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_from_base_for_u8() {
        assert_eq!(u8::from(Base::A), b'A');
        assert_eq!(u8::from(Base::C), b'C');
        assert_eq!(u8::from(Base::G), b'G');
        assert_eq!(u8::from(Base::T), b'T');
        assert_eq!(u8::from(Base::N), b'N');
    }
}
