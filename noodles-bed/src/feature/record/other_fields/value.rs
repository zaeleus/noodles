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

impl<'a> From<Value<'a>> for crate::feature::record_buf::other_fields::Value {
    fn from(value: Value<'a>) -> Self {
        match value {
            Value::Int64(n) => Self::Int64(n),
            Value::UInt64(n) => Self::UInt64(n),
            Value::Float64(n) => Self::Float64(n),
            Value::Character(b) => Self::Character(b),
            Value::String(s) => Self::String(s.into()),
        }
    }
}
