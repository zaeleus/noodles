#[derive(Clone, Copy, Debug, Eq, PartialEq)]
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

impl From<Base> for char {
    fn from(base: Base) -> char {
        match base {
            Base::A => 'A',
            Base::C => 'C',
            Base::G => 'G',
            Base::T => 'T',
            Base::N => 'N',
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
    fn test_from_base_for_char() {
        assert_eq!(char::from(Base::A), 'A');
        assert_eq!(char::from(Base::C), 'C');
        assert_eq!(char::from(Base::G), 'G');
        assert_eq!(char::from(Base::T), 'T');
        assert_eq!(char::from(Base::N), 'N');
    }
}
