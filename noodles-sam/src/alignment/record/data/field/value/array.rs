//! Alignment record data field array value.

mod subtype;
mod values;

use std::{fmt, io};

pub use self::{subtype::Subtype, values::Values};

/// An alignment record data field array value.
pub enum Array<'a> {
    /// An 8-bit integer array (`B:c`).
    Int8(Box<dyn Values<'a, i8> + 'a>),
    /// An 8-bit unsigned integer array (`B:C`).
    UInt8(Box<dyn Values<'a, u8> + 'a>),
    /// A 16-bit integer array (`B:s`).
    Int16(Box<dyn Values<'a, i16> + 'a>),
    /// A 16-bit unsigned integer array (`B:S`).
    UInt16(Box<dyn Values<'a, u16> + 'a>),
    /// A 32-bit integer array (`B:i`).
    Int32(Box<dyn Values<'a, i32> + 'a>),
    /// A 32-bit unsigned integer array (`B:I`).
    UInt32(Box<dyn Values<'a, u32> + 'a>),
    /// A single-precision floating-point array (`B:f`).
    Float(Box<dyn Values<'a, f32> + 'a>),
}

impl Array<'_> {
    /// Returns the array value type.
    pub fn subtype(&self) -> Subtype {
        match self {
            Array::Int8(_) => Subtype::Int8,
            Array::UInt8(_) => Subtype::UInt8,
            Array::Int16(_) => Subtype::Int16,
            Array::UInt16(_) => Subtype::UInt16,
            Array::Int32(_) => Subtype::Int32,
            Array::UInt32(_) => Subtype::UInt32,
            Array::Float(_) => Subtype::Float,
        }
    }
}

impl fmt::Debug for Array<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Int8(values) => f.debug_list().entries(values.iter()).finish(),
            Self::UInt8(values) => f.debug_list().entries(values.iter()).finish(),
            Self::Int16(values) => f.debug_list().entries(values.iter()).finish(),
            Self::UInt16(values) => f.debug_list().entries(values.iter()).finish(),
            Self::Int32(values) => f.debug_list().entries(values.iter()).finish(),
            Self::UInt32(values) => f.debug_list().entries(values.iter()).finish(),
            Self::Float(values) => f.debug_list().entries(values.iter()).finish(),
        }
    }
}

impl<'a> TryFrom<Array<'a>> for crate::alignment::record_buf::data::field::value::Array {
    type Error = io::Error;

    fn try_from(array: Array<'a>) -> Result<Self, Self::Error> {
        match array {
            Array::Int8(values) => values.iter().collect::<Result<_, _>>().map(Self::Int8),
            Array::UInt8(values) => values.iter().collect::<Result<_, _>>().map(Self::UInt8),
            Array::Int16(values) => values.iter().collect::<Result<_, _>>().map(Self::Int16),
            Array::UInt16(values) => values.iter().collect::<Result<_, _>>().map(Self::UInt16),
            Array::Int32(values) => values.iter().collect::<Result<_, _>>().map(Self::Int32),
            Array::UInt32(values) => values.iter().collect::<Result<_, _>>().map(Self::UInt32),
            Array::Float(values) => values.iter().collect::<Result<_, _>>().map(Self::Float),
        }
    }
}
