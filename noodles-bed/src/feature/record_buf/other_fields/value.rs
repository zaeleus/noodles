use bstr::BString;

/// A feature record other field value buffer.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    /// A 64-bit signed integer.
    Int64(i64),
    /// A 64-bit unsigned integer.
    UInt64(u64),
    /// A double-precision floating-point.
    Float64(f64),
    /// An ASCII character.
    Character(u8),
    /// A byte string.
    String(BString),
}

impl From<i64> for Value {
    fn from(n: i64) -> Self {
        Self::Int64(n)
    }
}

impl From<u64> for Value {
    fn from(n: u64) -> Self {
        Self::UInt64(n)
    }
}

impl From<f64> for Value {
    fn from(n: f64) -> Self {
        Self::Float64(n)
    }
}

impl From<&str> for Value {
    fn from(s: &str) -> Self {
        Self::from(String::from(s))
    }
}

impl From<String> for Value {
    fn from(s: String) -> Self {
        Self::String(s.into())
    }
}

impl<'a> From<&'a Value> for crate::feature::record::other_fields::Value<'a> {
    fn from(value_buf: &'a Value) -> Self {
        match value_buf {
            Value::Int64(n) => Self::Int64(*n),
            Value::UInt64(n) => Self::UInt64(*n),
            Value::Float64(n) => Self::Float64(*n),
            Value::Character(c) => Self::Character(*c),
            Value::String(s) => Self::String(s.as_ref()),
        }
    }
}
