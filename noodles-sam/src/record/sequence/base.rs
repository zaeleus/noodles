use std::{convert::TryFrom, error, fmt};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Base {
    A,
    B,
    C,
    D,
    E,
    F,
    G,
    H,
    I,
    J,
    K,
    L,
    M,
    N,
    O,
    P,
    Q,
    R,
    S,
    T,
    U,
    V,
    W,
    X,
    Y,
    Z,
}

impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", char::from(*self))
    }
}

#[derive(Debug)]
pub struct TryFromCharError(char);

impl error::Error for TryFromCharError {}

impl fmt::Display for TryFromCharError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid base: expected {{A..=Z}}, got {}", self.0)
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
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_char_for_base() -> Result<(), TryFromCharError> {
        assert_eq!(Base::try_from('A')?, Base::A);
        assert_eq!(Base::try_from('B')?, Base::B);
        assert_eq!(Base::try_from('C')?, Base::C);
        assert_eq!(Base::try_from('D')?, Base::D);
        assert_eq!(Base::try_from('E')?, Base::E);
        assert_eq!(Base::try_from('F')?, Base::F);
        assert_eq!(Base::try_from('G')?, Base::G);
        assert_eq!(Base::try_from('H')?, Base::H);
        assert_eq!(Base::try_from('I')?, Base::I);
        assert_eq!(Base::try_from('J')?, Base::J);
        assert_eq!(Base::try_from('K')?, Base::K);
        assert_eq!(Base::try_from('L')?, Base::L);
        assert_eq!(Base::try_from('M')?, Base::M);
        assert_eq!(Base::try_from('N')?, Base::N);
        assert_eq!(Base::try_from('O')?, Base::O);
        assert_eq!(Base::try_from('P')?, Base::P);
        assert_eq!(Base::try_from('Q')?, Base::Q);
        assert_eq!(Base::try_from('R')?, Base::R);
        assert_eq!(Base::try_from('S')?, Base::S);
        assert_eq!(Base::try_from('T')?, Base::T);
        assert_eq!(Base::try_from('U')?, Base::U);
        assert_eq!(Base::try_from('V')?, Base::V);
        assert_eq!(Base::try_from('W')?, Base::W);
        assert_eq!(Base::try_from('X')?, Base::X);
        assert_eq!(Base::try_from('Y')?, Base::Y);
        assert_eq!(Base::try_from('Z')?, Base::Z);

        assert!(Base::try_from('*').is_err());

        Ok(())
    }
}
