// § 6.3.3.2 "Type encoding: Integers" (2024-10-09)
const MISSING: i16 = i16::MIN;
const END_OF_VECTOR: i16 = i16::MIN + 1;
const RESERVED_0: i16 = i16::MIN + 2;
const RESERVED_5: i16 = i16::MIN + 7;

#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum Int16 {
    Value(i16),
    Missing,
    EndOfVector,
    Reserved(i16),
}

impl Int16 {
    pub const MIN_VALUE: i16 = i16::MIN + 8;
    pub const MAX_VALUE: i16 = i16::MAX;
}

impl From<i16> for Int16 {
    fn from(value: i16) -> Self {
        match value {
            MISSING => Self::Missing,
            END_OF_VECTOR => Self::EndOfVector,
            RESERVED_0..=RESERVED_5 => Self::Reserved(value),
            _ => Self::Value(value),
        }
    }
}

impl From<Int16> for i16 {
    fn from(value: Int16) -> Self {
        match value {
            Int16::Missing => MISSING,
            Int16::EndOfVector => END_OF_VECTOR,
            Int16::Value(n) | Int16::Reserved(n) => n,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_i16_for_int16() {
        assert_eq!(Int16::from(-32768), Int16::Missing);
        assert_eq!(Int16::from(-32767), Int16::EndOfVector);
        assert_eq!(Int16::from(-32766), Int16::Reserved(-32766));
        assert_eq!(Int16::from(-32765), Int16::Reserved(-32765));
        assert_eq!(Int16::from(-32764), Int16::Reserved(-32764));
        assert_eq!(Int16::from(-32763), Int16::Reserved(-32763));
        assert_eq!(Int16::from(-32762), Int16::Reserved(-32762));
        assert_eq!(Int16::from(-32761), Int16::Reserved(-32761));
        assert_eq!(Int16::from(-32760), Int16::Value(-32760));
        assert_eq!(Int16::from(0), Int16::Value(0));
        assert_eq!(Int16::from(i16::MAX), Int16::Value(i16::MAX));
    }

    #[test]
    fn test_from_int16_for_i16() {
        assert_eq!(i16::from(Int16::Missing), -32768);
        assert_eq!(i16::from(Int16::EndOfVector), -32767);
        assert_eq!(i16::from(Int16::Reserved(-32766)), -32766);
        assert_eq!(i16::from(Int16::Reserved(-32765)), -32765);
        assert_eq!(i16::from(Int16::Reserved(-32764)), -32764);
        assert_eq!(i16::from(Int16::Reserved(-32763)), -32763);
        assert_eq!(i16::from(Int16::Reserved(-32762)), -32762);
        assert_eq!(i16::from(Int16::Reserved(-32761)), -32761);
        assert_eq!(i16::from(Int16::Value(-32760)), -32760);
        assert_eq!(i16::from(Int16::Value(0)), 0);
        assert_eq!(i16::from(Int16::Value(i16::MAX)), i16::MAX);
    }
}
