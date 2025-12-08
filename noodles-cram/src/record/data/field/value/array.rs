mod values;

use noodles_sam as sam;

pub use self::values::Values;

#[derive(Clone, Debug, PartialEq)]
pub enum Array<'c> {
    Int8(Values<'c, i8>),
    UInt8(Values<'c, u8>),
    Int16(Values<'c, i16>),
    UInt16(Values<'c, u16>),
    Int32(Values<'c, i32>),
    UInt32(Values<'c, u32>),
    Float(Values<'c, f32>),
}

impl<'c> From<Array<'c>> for sam::alignment::record::data::field::value::Array<'c> {
    fn from(array: Array<'c>) -> Self {
        match array {
            Array::Int8(values) => Self::Int8(Box::new(values)),
            Array::UInt8(values) => Self::UInt8(Box::new(values)),
            Array::Int16(values) => Self::Int16(Box::new(values)),
            Array::UInt16(values) => Self::UInt16(Box::new(values)),
            Array::Int32(values) => Self::Int32(Box::new(values)),
            Array::UInt32(values) => Self::UInt32(Box::new(values)),
            Array::Float(values) => Self::Float(Box::new(values)),
        }
    }
}
