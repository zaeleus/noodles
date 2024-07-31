use bstr::BStr;

/// A feature record other field value.
#[derive(Debug)]
pub enum Value<'a> {
    /// A 64-bit signed integer.
    Int64(i64),
    /// A 64-bit unsigned integer.
    UInt64(u64),
    /// A double-precision floating-point.
    Float64(f64),
    /// An ASCII character.
    Character(u8),
    /// A byte string.
    String(&'a BStr),
}
