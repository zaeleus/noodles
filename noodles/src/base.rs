#[derive(Debug, Eq, PartialEq)]
pub enum Base {
    A,
    C,
    G,
    T,
    U,
}

impl Base {
    pub fn from_char(c: char) -> Result<Base, ()> {
        match c {
            'A' => Ok(Base::A),
            'C' => Ok(Base::C),
            'G' => Ok(Base::G),
            'T' => Ok(Base::T),
            'U' => Ok(Base::U),
            _ => Err(()),
        }
    }

    pub fn complement(&self) -> Base {
        match self {
            Base::A => Base::T,
            Base::C => Base::G,
            Base::G => Base::C,
            Base::T => Base::A,
            Base::U => Base::A,
        }
    }

    pub fn description(&self) -> &'static str {
        match self {
            Base::A => "adenine",
            Base::C => "cytosine",
            Base::G => "guanine",
            Base::T => "thymine",
            Base::U => "uracil",
        }
    }

    pub fn symbol(&self) -> char {
        match self {
            Base::A => 'A',
            Base::C => 'C',
            Base::G => 'G',
            Base::T => 'T',
            Base::U => 'U',
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
