pub mod array;

use noodles_sam as sam;

pub use self::array::Array;
use bstr::BStr;

#[derive(Clone, Debug, PartialEq)]
pub enum Value<'c> {
    Character(u8),
    Int8(i8),
    UInt8(u8),
    Int16(i16),
    UInt16(u16),
    Int32(i32),
    UInt32(u32),
    Float(f32),
    String(&'c BStr),
    Hex(&'c BStr),
    Array(Array<'c>),
}

impl<'r: 'c, 'c: 'r> From<&'r Value<'c>> for sam::alignment::record::data::field::Value<'c> {
    fn from(value: &'r Value<'c>) -> Self {
        match value {
            Value::Character(c) => Self::Character(*c),
            Value::Int8(n) => Self::Int8(*n),
            Value::UInt8(n) => Self::UInt8(*n),
            Value::Int16(n) => Self::Int16(*n),
            Value::UInt16(n) => Self::UInt16(*n),
            Value::Int32(n) => Self::Int32(*n),
            Value::UInt32(n) => Self::UInt32(*n),
            Value::Float(n) => Self::Float(*n),
            Value::String(s) => Self::String(s),
            Value::Hex(s) => Self::Hex(s),
            Value::Array(array) => Self::Array(array.into()),
        }
    }
}
