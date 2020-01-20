#[derive(Debug, Eq, PartialEq)]
pub enum Base {
    A,
    C,
    G,
    T,
    U,
}

impl Base {
    pub fn from_char(c: char) -> Result<Self, ()> {
        match c {
            'A' => Ok(Self::A),
            'C' => Ok(Self::C),
            'G' => Ok(Self::G),
            'T' => Ok(Self::T),
            'U' => Ok(Self::U),
            _ => Err(()),
        }
    }

    pub fn complement(&self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
            Self::U => Self::A,
        }
    }

    pub fn description(&self) -> &'static str {
        match self {
            Self::A => "adenine",
            Self::C => "cytosine",
            Self::G => "guanine",
            Self::T => "thymine",
            Self::U => "uracil",
        }
    }

    pub fn symbol(&self) -> char {
        match self {
            Self::A => 'A',
            Self::C => 'C',
            Self::G => 'G',
            Self::T => 'T',
            Self::U => 'U',
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Base;

    #[test]
    fn test_from_char() {
        assert_eq!(Base::from_char('A'), Ok(Base::A));
        assert_eq!(Base::from_char('C'), Ok(Base::C));
        assert_eq!(Base::from_char('G'), Ok(Base::G));
        assert_eq!(Base::from_char('T'), Ok(Base::T));
        assert_eq!(Base::from_char('U'), Ok(Base::U));

        assert!(Base::from_char('Q').is_err());
    }

    #[test]
    fn test_complement() {
        assert_eq!(Base::A.complement(), Base::T);
        assert_eq!(Base::C.complement(), Base::G);
        assert_eq!(Base::G.complement(), Base::C);
        assert_eq!(Base::T.complement(), Base::A);
        assert_eq!(Base::U.complement(), Base::A);
    }

    #[test]
    fn test_description() {
        assert_eq!(Base::A.description(), "adenine");
        assert_eq!(Base::C.description(), "cytosine");
        assert_eq!(Base::G.description(), "guanine");
        assert_eq!(Base::T.description(), "thymine");
        assert_eq!(Base::U.description(), "uracil");
    }

    #[test]
    fn test_symbol() {
        assert_eq!(Base::A.symbol(), 'A');
        assert_eq!(Base::C.symbol(), 'C');
        assert_eq!(Base::G.symbol(), 'G');
        assert_eq!(Base::T.symbol(), 'T');
        assert_eq!(Base::U.symbol(), 'U');
    }
}
