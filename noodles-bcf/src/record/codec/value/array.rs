#[derive(Clone, Debug, PartialEq)]
pub enum Array {
    Int8(Vec<i8>),
    Int16(Vec<i16>),
    Int32(Vec<i32>),
    Float(Vec<f32>),
}
