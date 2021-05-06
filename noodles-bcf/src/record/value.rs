mod float;
mod int16;
mod int32;
mod int8;
mod ty;

pub use self::{float::Float, int16::Int16, int32::Int32, int8::Int8, ty::Type};

#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    Int8(Option<Int8>),
    Int8Array(Vec<i8>),
    Int16(Option<Int16>),
    Int16Array(Vec<i16>),
    Int32(Option<Int32>),
    Int32Array(Vec<i32>),
    Float(Option<Float>),
    FloatArray(Vec<f32>),
    String(Option<String>),
}
