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
        match n.to_ascii_uppercase() {
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

impl From<Base> for u8 {
    fn from(base: Base) -> Self {
        match base {
            Base::A => b'A',
            Base::B => b'B',
            Base::C => b'C',
            Base::D => b'D',
            Base::E => b'E',
            Base::F => b'F',
            Base::G => b'G',
            Base::H => b'H',
            Base::I => b'I',
            Base::J => b'J',
            Base::K => b'K',
            Base::L => b'L',
            Base::M => b'M',
            Base::N => b'N',
            Base::O => b'O',
            Base::P => b'P',
            Base::Q => b'Q',
            Base::R => b'R',
            Base::S => b'S',
            Base::T => b'T',
            Base::U => b'U',
            Base::V => b'V',
            Base::W => b'W',
            Base::X => b'X',
            Base::Y => b'Y',
            Base::Z => b'Z',
            Base::Eq => b'=',
        }
    }
}

#[cfg(test)]
mod tests {
    #[rustfmt::skip]
    static ALPHA_BASES: &[Base; 26] = &[
        Base::A, Base::B, Base::C, Base::D, Base::E, Base::F, Base::G, Base::H, Base::I, Base::J,
        Base::K, Base::L, Base::M, Base::N, Base::O, Base::P, Base::Q, Base::R, Base::S, Base::T,
        Base::U, Base::V, Base::W, Base::X, Base::Y, Base::Z,
    ];

    use super::*;

    #[test]
    fn test_try_from_char_for_base() {
        for (c, &expected) in ('A'..='Z').zip(ALPHA_BASES) {
            assert_eq!(Base::try_from(c), Ok(expected));
            assert_eq!(Base::try_from(c.to_ascii_lowercase()), Ok(expected));
        }

        assert_eq!(Base::try_from('='), Ok(Base::Eq));

        assert_eq!(Base::try_from('*'), Err(TryFromCharError('*')));
    }

    #[test]
    fn test_try_from_u8_for_base() {
        for (c, &expected) in (b'A'..=b'Z').zip(ALPHA_BASES) {
            assert_eq!(Base::try_from(c), Ok(expected));
            assert_eq!(Base::try_from(c.to_ascii_lowercase()), Ok(expected));
        }

        assert_eq!(Base::try_from(b'='), Ok(Base::Eq));

        assert_eq!(Base::try_from(b'*'), Err(TryFromCharError('*')));
    }

    #[test]
    fn test_from_base_for_u8() {
        for (&base, expected) in ALPHA_BASES.iter().zip(b'A'..=b'Z') {
            assert_eq!(u8::from(base), expected);
        }

        assert_eq!(u8::from(Base::Eq), b'=');
    }
}
