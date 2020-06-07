use std::ops::Deref;

const UNMAPPED: i32 = 0;
const MIN: i32 = 1;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Position(Option<i32>);

impl From<i32> for Position {
    fn from(n: i32) -> Self {
        if n < MIN {
            Self(None)
        } else {
            Self(Some(n))
        }
    }
}

impl From<Position> for i32 {
    fn from(position: Position) -> Self {
        match *position {
            Some(n) => n,
            None => UNMAPPED,
        }
    }
}

impl Deref for Position {
    type Target = Option<i32>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_i32_for_position() {
        assert_eq!(*Position::from(0), None);
        assert_eq!(*Position::from(13), Some(13));
    }

    #[test]
    fn test_from_position_for_i32() {
        assert_eq!(i32::from(Position::from(0)), 0);
        assert_eq!(i32::from(Position::from(13)), 13);
    }
}
