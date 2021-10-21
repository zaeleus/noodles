use std::{error, fmt};

#[derive(Debug, Eq, PartialEq)]
pub struct TryFromByteError(u8);

impl fmt::Display for TryFromByteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid rANS order: expected 0 or 1, got {}", self.0)
    }
}

impl error::Error for TryFromByteError {}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Order {
    Zero,
    One,
}

impl TryFrom<u8> for Order {
    type Error = TryFromByteError;
    fn try_from(b: u8) -> Result<Self, Self::Error> {
        match b {
            0 => Ok(Self::Zero),
            1 => Ok(Self::One),
            _ => Err(TryFromByteError(b)),
        }
    }
}

impl From<Order> for u8 {
    fn from(order: Order) -> Self {
        order as Self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_u8_for_order() {
        assert_eq!(Order::try_from(0), Ok(Order::Zero));
        assert_eq!(Order::try_from(1), Ok(Order::One));
        assert_eq!(Order::try_from(2), Err(TryFromByteError(2)));
    }

    #[test]
    fn test_from_order_for_u8() {
        assert_eq!(u8::from(Order::Zero), 0);
        assert_eq!(u8::from(Order::One), 1);
    }
}
