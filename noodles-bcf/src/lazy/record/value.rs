mod array;
mod float;
mod int16;
mod int32;
mod int8;
mod ty;

pub use self::{array::Array, float::Float, int16::Int16, int32::Int32, int8::Int8, ty::Type};

#[derive(Clone, Debug, PartialEq)]
pub enum Value<'a> {
    Int8(Option<Int8>),
    Int16(Option<Int16>),
    Int32(Option<Int32>),
    Float(Option<Float>),
    String(Option<&'a str>),
    Array(Array),
}
