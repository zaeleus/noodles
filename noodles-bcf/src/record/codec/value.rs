pub mod array;
mod float;
mod int16;
mod int32;
mod int8;
mod ty;

pub use self::{array::Array, float::Float, int16::Int16, int32::Int32, int8::Int8, ty::Type};

#[derive(Debug)]
pub enum Value<'a> {
    Int8(Option<Int8>),
    Int16(Option<Int16>),
    Int32(Option<Int32>),
    Float(Option<Float>),
    String(Option<&'a str>),
    Array(Array<'a>),
}

impl<'a> Value<'a> {
    pub fn as_int(&self) -> Option<i32> {
        match self {
            Self::Int8(Some(Int8::Value(n))) => Some(i32::from(*n)),
            Self::Int16(Some(Int16::Value(n))) => Some(i32::from(*n)),
            Self::Int32(Some(Int32::Value(n))) => Some(*n),
            _ => None,
        }
    }
}
