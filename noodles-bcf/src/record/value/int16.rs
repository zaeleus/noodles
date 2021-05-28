#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum Int16 {
    Value(i16),
    Missing,
    EndOfVector,
    Reserved(i16),
}

impl Int16 {
    /// The smallest value that can be represented by [`Self::Value`].
    pub const MIN_VALUE: i16 = i16::MIN + 8;

    /// The largest value that can be represented by [`Self::Value`].
    pub const MAX_VALUE: i16 = i16::MAX;
}

impl From<i16> for Int16 {
    fn from(value: i16) -> Self {
        match value as u16 {
            0x8000 => Self::Missing,
            0x8001 => Self::EndOfVector,
            0x8002..=0x8007 => Self::Reserved(value),
            _ => Self::Value(value),
        }
    }
}

impl From<Int16> for i16 {
    fn from(value: Int16) -> Self {
        match value {
            Int16::Missing => -32768,
            Int16::EndOfVector => -32767,
            Int16::Value(n) | Int16::Reserved(n) => n,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_i16_for_int16() {
        assert_eq!(Int16::from(0), Int16::Value(0));
        assert_eq!(Int16::from(-32768), Int16::Missing);
        assert_eq!(Int16::from(-32767), Int16::EndOfVector);
        assert_eq!(Int16::from(-32766), Int16::Reserved(-32766));
        assert_eq!(Int16::from(-32765), Int16::Reserved(-32765));
        assert_eq!(Int16::from(-32764), Int16::Reserved(-32764));
        assert_eq!(Int16::from(-32763), Int16::Reserved(-32763));
        assert_eq!(Int16::from(-32762), Int16::Reserved(-32762));
        assert_eq!(Int16::from(-32761), Int16::Reserved(-32761));
    }

    #[test]
    fn test_from_int16_for_i16() {
        assert_eq!(i16::from(Int16::Value(0)), 0);
        assert_eq!(i16::from(Int16::Missing), -32768);
        assert_eq!(i16::from(Int16::EndOfVector), -32767);
        assert_eq!(i16::from(Int16::Reserved(-32766)), -32766);
        assert_eq!(i16::from(Int16::Reserved(-32765)), -32765);
        assert_eq!(i16::from(Int16::Reserved(-32764)), -32764);
        assert_eq!(i16::from(Int16::Reserved(-32763)), -32763);
        assert_eq!(i16::from(Int16::Reserved(-32762)), -32762);
        assert_eq!(i16::from(Int16::Reserved(-32761)), -32761);
    }
}
