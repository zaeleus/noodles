use std::ops::Deref;

const UNMAPPED_POSITION: u32 = 0;

#[derive(Clone, Copy, Debug, Default)]
pub struct Position(Option<u32>);

impl From<u32> for Position {
    fn from(n: u32) -> Self {
        if n == UNMAPPED_POSITION {
            Self(None)
        } else {
            Self(Some(n))
        }
    }
}

impl From<Position> for u32 {
    fn from(position: Position) -> Self {
        match *position {
            Some(n) => n,
            None => UNMAPPED_POSITION,
        }
    }
}

impl Deref for Position {
    type Target = Option<u32>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_u32_for_position() {
        assert_eq!(*Position::from(0), None);
        assert_eq!(*Position::from(13), Some(13));
    }

    #[test]
    fn test_from_position_for_u32() {
        assert_eq!(u32::from(Position::from(0)), 0);
        assert_eq!(u32::from(Position::from(13)), 13);
    }
}
