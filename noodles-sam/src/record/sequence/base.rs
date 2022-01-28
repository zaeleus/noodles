use std::{
    error,
    fmt::{self, Write},
};

/// A SAM record sequence base.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Base {
    /// Adenine.
    A,
    /// Not A.
    B,
    /// Cytosine.
    C,
    /// Not C.
    D,
    /// Undefined (`E`).
    E,
    /// Undefined (`F`).
    F,
    /// Guanine.
    G,
    /// Not G.
    H,
    /// Undefined (`I`).
    I,
    /// Undefined (`J`).
    J,
    /// Keto.
    K,
    /// Undefined (`L`).
    L,
    /// Amino.
    M,
    /// Any base.
    N,
    /// Undefined (`O`).
    O,
    /// Purine.
    P,
    /// Undefined (`Q`).
    Q,
    /// Undefined (`R`).
    R,
    /// Strong.
    S,
    /// Thymine.
    T,
    /// Uracil.
    U,
    /// Not T.
    V,
    /// Weak.
    W,
    /// Undefined (`X`).
    X,
    /// Pyrimidine.
    Y,
    /// Zero.
    Z,
    /// Equal to the reference base.
    Eq,
}

impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_char(char::from(*self))
    }
}

/// An error returned when the conversion from a character to a SAM record sequence base fails.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromCharError(char);

impl error::Error for TryFromCharError {}

impl fmt::Display for TryFromCharError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "expected {{A..=Z, =}}, got {}", self.0)
    }
}

impl TryFrom<char> for Base {
    type Error = TryFromCharError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            'A' => Ok(Self::A),
            'B' => Ok(Self::B),
            'C' => Ok(Self::C),
            'D' => Ok(Self::D),
            'E' => Ok(Self::E),
            'F' => Ok(Self::F),
            'G' => Ok(Self::G),
            'H' => Ok(Self::H),
            'I' => Ok(Self::I),
            'J' => Ok(Self::J),
            'K' => Ok(Self::K),
            'L' => Ok(Self::L),
            'M' => Ok(Self::M),
            'N' => Ok(Self::N),
            'O' => Ok(Self::O),
            'P' => Ok(Self::P),
            'Q' => Ok(Self::Q),
            'R' => Ok(Self::R),
            'S' => Ok(Self::S),
            'T' => Ok(Self::T),
            'U' => Ok(Self::U),
            'V' => Ok(Self::V),
            'W' => Ok(Self::W),
            'X' => Ok(Self::X),
            'Y' => Ok(Self::Y),
            'Z' => Ok(Self::Z),
            '=' => Ok(Self::Eq),
            _ => Err(TryFromCharError(c)),
        }
    }
}

impl From<Base> for char {
    fn from(base: Base) -> Self {
        match base {
            Base::A => 'A',
            Base::B => 'B',
            Base::C => 'C',
            Base::D => 'D',
            Base::E => 'E',
            Base::F => 'F',
            Base::G => 'G',
            Base::H => 'H',
            Base::I => 'I',
            Base::J => 'J',
            Base::K => 'K',
            Base::L => 'L',
            Base::M => 'M',
            Base::N => 'N',
            Base::O => 'O',
            Base::P => 'P',
            Base::Q => 'Q',
            Base::R => 'R',
            Base::S => 'S',
            Base::T => 'T',
            Base::U => 'U',
            Base::V => 'V',
            Base::W => 'W',
            Base::X => 'X',
            Base::Y => 'Y',
            Base::Z => 'Z',
            Base::Eq => '=',
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_char_for_base() {
        assert_eq!(Base::try_from('A'), Ok(Base::A));
        assert_eq!(Base::try_from('B'), Ok(Base::B));
        assert_eq!(Base::try_from('C'), Ok(Base::C));
        assert_eq!(Base::try_from('D'), Ok(Base::D));
        assert_eq!(Base::try_from('E'), Ok(Base::E));
        assert_eq!(Base::try_from('F'), Ok(Base::F));
        assert_eq!(Base::try_from('G'), Ok(Base::G));
        assert_eq!(Base::try_from('H'), Ok(Base::H));
        assert_eq!(Base::try_from('I'), Ok(Base::I));
        assert_eq!(Base::try_from('J'), Ok(Base::J));
        assert_eq!(Base::try_from('K'), Ok(Base::K));
        assert_eq!(Base::try_from('L'), Ok(Base::L));
        assert_eq!(Base::try_from('M'), Ok(Base::M));
        assert_eq!(Base::try_from('N'), Ok(Base::N));
        assert_eq!(Base::try_from('O'), Ok(Base::O));
        assert_eq!(Base::try_from('P'), Ok(Base::P));
        assert_eq!(Base::try_from('Q'), Ok(Base::Q));
        assert_eq!(Base::try_from('R'), Ok(Base::R));
        assert_eq!(Base::try_from('S'), Ok(Base::S));
        assert_eq!(Base::try_from('T'), Ok(Base::T));
        assert_eq!(Base::try_from('U'), Ok(Base::U));
        assert_eq!(Base::try_from('V'), Ok(Base::V));
        assert_eq!(Base::try_from('W'), Ok(Base::W));
        assert_eq!(Base::try_from('X'), Ok(Base::X));
        assert_eq!(Base::try_from('Y'), Ok(Base::Y));
        assert_eq!(Base::try_from('Z'), Ok(Base::Z));
        assert_eq!(Base::try_from('='), Ok(Base::Eq));

        assert_eq!(Base::try_from('*'), Err(TryFromCharError('*')));
    }
}
