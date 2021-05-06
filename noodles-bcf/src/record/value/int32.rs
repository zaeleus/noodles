#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum Int32 {
    Value(i32),
    Missing,
    EndOfVector,
    Reserved(i32),
}

impl From<i32> for Int32 {
    fn from(value: i32) -> Self {
        match value as u32 {
            0x80000000 => Self::Missing,
            0x80000001 => Self::EndOfVector,
            0x80000002..=0x80000007 => Self::Reserved(value),
            _ => Self::Value(value),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_i32_for_int32() {
        assert_eq!(Int32::from(0), Int32::Value(0));
        assert_eq!(Int32::from(-2147483648), Int32::Missing);
        assert_eq!(Int32::from(-2147483647), Int32::EndOfVector);
        assert_eq!(Int32::from(-2147483646), Int32::Reserved(-2147483646));
        assert_eq!(Int32::from(-2147483645), Int32::Reserved(-2147483645));
        assert_eq!(Int32::from(-2147483644), Int32::Reserved(-2147483644));
        assert_eq!(Int32::from(-2147483643), Int32::Reserved(-2147483643));
        assert_eq!(Int32::from(-2147483642), Int32::Reserved(-2147483642));
        assert_eq!(Int32::from(-2147483641), Int32::Reserved(-2147483641));
    }
}
