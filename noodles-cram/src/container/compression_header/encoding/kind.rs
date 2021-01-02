#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    Null,
    External,
    Golomb,
    Huffman,
    ByteArrayLen,
    ByteArrayStop,
    Beta,
    Subexp,
    GolombRice,
    Gamma,
}

impl From<Kind> for i32 {
    fn from(kind: Kind) -> Self {
        kind as Self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_kind_for_i32() {
        assert_eq!(i32::from(Kind::Null), 0);
        assert_eq!(i32::from(Kind::External), 1);
        assert_eq!(i32::from(Kind::Golomb), 2);
        assert_eq!(i32::from(Kind::Huffman), 3);
        assert_eq!(i32::from(Kind::ByteArrayLen), 4);
        assert_eq!(i32::from(Kind::ByteArrayStop), 5);
        assert_eq!(i32::from(Kind::Beta), 6);
        assert_eq!(i32::from(Kind::Subexp), 7);
        assert_eq!(i32::from(Kind::GolombRice), 8);
        assert_eq!(i32::from(Kind::Gamma), 9);
    }
}
