//! Alignment record data field array value buffer.

use std::io;

use crate::alignment::record::data::field::value::array::Subtype;

/// An alignment record data field array value buffer.
#[derive(Clone, Debug, PartialEq)]
pub enum Array {
    /// An 8-bit integer array (`B:c`).
    Int8(Vec<i8>),
    /// An 8-bit unsigned integer array (`B:C`).
    UInt8(Vec<u8>),
    /// A 16-bit integer array (`B:s`).
    Int16(Vec<i16>),
    /// A 16-bit unsigned integer array (`B:S`).
    UInt16(Vec<u16>),
    /// A 32-bit integer array (`B:i`).
    Int32(Vec<i32>),
    /// A 32-bit unsigned integer array (`B:I`).
    UInt32(Vec<u32>),
    /// A single-precision floating-point array (`B:f`).
    Float(Vec<f32>),
}

impl Array {
    /// Returns the type of the array values.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::{
    ///     record::data::field::value::array::Subtype, record_buf::data::field::value::Array,
    /// };
    ///
    /// assert_eq!(Array::UInt8(Vec::new()).subtype(), Subtype::UInt8);
    /// ```
    pub fn subtype(&self) -> Subtype {
        match self {
            Self::Int8(_) => Subtype::Int8,
            Self::UInt8(_) => Subtype::UInt8,
            Self::Int16(_) => Subtype::Int16,
            Self::UInt16(_) => Subtype::UInt16,
            Self::Int32(_) => Subtype::Int32,
            Self::UInt32(_) => Subtype::UInt32,
            Self::Float(_) => Subtype::Float,
        }
    }
}

impl<'a> From<&'a Array> for crate::alignment::record::data::field::value::Array<'a> {
    fn from(array_buf: &'a Array) -> Self {
        match array_buf {
            Array::Int8(values) => Self::Int8(Box::new(Values::new(values))),
            Array::UInt8(values) => Self::UInt8(Box::new(Values::new(values))),
            Array::Int16(values) => Self::Int16(Box::new(Values::new(values))),
            Array::UInt16(values) => Self::UInt16(Box::new(Values::new(values))),
            Array::Int32(values) => Self::Int32(Box::new(Values::new(values))),
            Array::UInt32(values) => Self::UInt32(Box::new(Values::new(values))),
            Array::Float(values) => Self::Float(Box::new(Values::new(values))),
        }
    }
}

struct Values<'a, N>(&'a [N]);

impl<'a, N> Values<'a, N> {
    fn new(values: &'a [N]) -> Self {
        Self(values)
    }
}

impl<'a, N> crate::alignment::record::data::field::value::array::Values<'a, N> for Values<'a, N>
where
    N: Copy,
{
    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<N>> + '_> {
        Box::new(self.0.iter().copied().map(Ok))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_subtype() {
        assert_eq!(Array::Int8(Vec::new()).subtype(), Subtype::Int8);
        assert_eq!(Array::UInt8(Vec::new()).subtype(), Subtype::UInt8);
        assert_eq!(Array::Int16(Vec::new()).subtype(), Subtype::Int16);
        assert_eq!(Array::UInt16(Vec::new()).subtype(), Subtype::UInt16);
        assert_eq!(Array::Int32(Vec::new()).subtype(), Subtype::Int32);
        assert_eq!(Array::UInt32(Vec::new()).subtype(), Subtype::UInt32);
        assert_eq!(Array::Float(Vec::new()).subtype(), Subtype::Float);
    }
}
