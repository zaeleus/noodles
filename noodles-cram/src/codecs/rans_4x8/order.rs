/// The number of bytes of context used to compute frequencies.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Order {
    /// Order-0.
    Zero,
    /// Order-1.
    One,
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
    fn test_from_order_for_u8() {
        assert_eq!(u8::from(Order::Zero), 0);
        assert_eq!(u8::from(Order::One), 1);
    }
}
