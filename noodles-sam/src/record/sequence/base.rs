//! SAM record sequence base.

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
        u8::try_from(c)
            .map_err(|_| TryFromCharError(c))
            .and_then(Self::try_from)
    }
}

impl TryFrom<u8> for Base {
    type Error = TryFromCharError;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        match n {
            b'A' => Ok(Self::A),
            b'B' => Ok(Self::B),
            b'C' => Ok(Self::C),
            b'D' => Ok(Self::D),
            b'E' => Ok(Self::E),
            b'F' => Ok(Self::F),
            b'G' => Ok(Self::G),
            b'H' => Ok(Self::H),
            b'I' => Ok(Self::I),
            b'J' => Ok(Self::J),
            b'K' => Ok(Self::K),
            b'L' => Ok(Self::L),
            b'M' => Ok(Self::M),
            b'N' => Ok(Self::N),
            b'O' => Ok(Self::O),
            b'P' => Ok(Self::P),
            b'Q' => Ok(Self::Q),
            b'R' => Ok(Self::R),
            b'S' => Ok(Self::S),
            b'T' => Ok(Self::T),
            b'U' => Ok(Self::U),
            b'V' => Ok(Self::V),
            b'W' => Ok(Self::W),
            b'X' => Ok(Self::X),
            b'Y' => Ok(Self::Y),
            b'Z' => Ok(Self::Z),
            b'=' => Ok(Self::Eq),
            _ => Err(TryFromCharError(char::from(n))),
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

    #[test]
    fn test_try_from_u8_for_base() {
        assert_eq!(Base::try_from(b'A'), Ok(Base::A));
        assert_eq!(Base::try_from(b'B'), Ok(Base::B));
        assert_eq!(Base::try_from(b'C'), Ok(Base::C));
        assert_eq!(Base::try_from(b'D'), Ok(Base::D));
        assert_eq!(Base::try_from(b'E'), Ok(Base::E));
        assert_eq!(Base::try_from(b'F'), Ok(Base::F));
        assert_eq!(Base::try_from(b'G'), Ok(Base::G));
        assert_eq!(Base::try_from(b'H'), Ok(Base::H));
        assert_eq!(Base::try_from(b'I'), Ok(Base::I));
        assert_eq!(Base::try_from(b'J'), Ok(Base::J));
        assert_eq!(Base::try_from(b'K'), Ok(Base::K));
        assert_eq!(Base::try_from(b'L'), Ok(Base::L));
        assert_eq!(Base::try_from(b'M'), Ok(Base::M));
        assert_eq!(Base::try_from(b'N'), Ok(Base::N));
        assert_eq!(Base::try_from(b'O'), Ok(Base::O));
        assert_eq!(Base::try_from(b'P'), Ok(Base::P));
        assert_eq!(Base::try_from(b'Q'), Ok(Base::Q));
        assert_eq!(Base::try_from(b'R'), Ok(Base::R));
        assert_eq!(Base::try_from(b'S'), Ok(Base::S));
        assert_eq!(Base::try_from(b'T'), Ok(Base::T));
        assert_eq!(Base::try_from(b'U'), Ok(Base::U));
        assert_eq!(Base::try_from(b'V'), Ok(Base::V));
        assert_eq!(Base::try_from(b'W'), Ok(Base::W));
        assert_eq!(Base::try_from(b'X'), Ok(Base::X));
        assert_eq!(Base::try_from(b'Y'), Ok(Base::Y));
        assert_eq!(Base::try_from(b'Z'), Ok(Base::Z));
        assert_eq!(Base::try_from(b'='), Ok(Base::Eq));

        assert_eq!(Base::try_from(b'*'), Err(TryFromCharError('*')));
    }
}
