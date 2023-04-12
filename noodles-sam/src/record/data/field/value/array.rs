//! SAM record data field array.

pub mod subtype;

pub use self::subtype::Subtype;

use std::fmt;

/// A SAM record data field array value.
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
    /// Returns the type of the array.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::{value::{Array, array::Subtype}, Value};
    /// assert_eq!(Array::UInt8(vec![0]).subtype(), Subtype::UInt8);
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

impl fmt::Display for Array {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Int8(values) => {
                char::from(Subtype::Int8).fmt(f)?;

                for value in values {
                    write!(f, ",{value}")?;
                }

                Ok(())
            }
            Self::UInt8(values) => {
                char::from(Subtype::UInt8).fmt(f)?;

                for value in values {
                    write!(f, ",{value}")?;
                }

                Ok(())
            }
            Self::Int16(values) => {
                char::from(Subtype::Int16).fmt(f)?;

                for value in values {
                    write!(f, ",{value}")?;
                }

                Ok(())
            }
            Self::UInt16(values) => {
                char::from(Subtype::UInt16).fmt(f)?;

                for value in values {
                    write!(f, ",{value}")?;
                }

                Ok(())
            }
            Self::Int32(values) => {
                char::from(Subtype::Int32).fmt(f)?;

                for value in values {
                    write!(f, ",{value}")?;
                }

                Ok(())
            }
            Self::UInt32(values) => {
                char::from(Subtype::UInt32).fmt(f)?;

                for value in values {
                    write!(f, ",{value}")?;
                }

                Ok(())
            }
            Self::Float(values) => {
                char::from(Subtype::Float).fmt(f)?;

                for value in values {
                    write!(f, ",{value}")?;
                }

                Ok(())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_subtype() {
        assert_eq!(Array::Int8(vec![0]).subtype(), Subtype::Int8);
        assert_eq!(Array::UInt8(vec![0]).subtype(), Subtype::UInt8);
        assert_eq!(Array::Int16(vec![0]).subtype(), Subtype::Int16);
        assert_eq!(Array::UInt16(vec![0]).subtype(), Subtype::UInt16);
        assert_eq!(Array::Int32(vec![0]).subtype(), Subtype::Int32);
        assert_eq!(Array::UInt32(vec![0]).subtype(), Subtype::UInt32);
        assert_eq!(Array::Float(vec![0.0]).subtype(), Subtype::Float);
    }

    #[test]
    fn test_fmt() {
        assert_eq!(Array::Int8(Vec::new()).to_string(), "c");
        assert_eq!(Array::Int8(vec![1, -2]).to_string(), "c,1,-2");

        assert_eq!(Array::UInt8(Vec::new()).to_string(), "C");
        assert_eq!(Array::UInt8(vec![3, 5]).to_string(), "C,3,5");

        assert_eq!(Array::Int16(Vec::new()).to_string(), "s");
        assert_eq!(Array::Int16(vec![8, -13]).to_string(), "s,8,-13");

        assert_eq!(Array::UInt16(Vec::new()).to_string(), "S");
        assert_eq!(Array::UInt16(vec![21, 34]).to_string(), "S,21,34");

        assert_eq!(Array::Int32(Vec::new()).to_string(), "i");
        assert_eq!(Array::Int32(vec![55, -89]).to_string(), "i,55,-89");

        assert_eq!(Array::UInt32(Vec::new()).to_string(), "I");
        assert_eq!(Array::UInt32(vec![144, 233]).to_string(), "I,144,233");

        assert_eq!(Array::Float(Vec::new()).to_string(), "f");
        assert_eq!(Array::Float(vec![0.0, 1.0]).to_string(), "f,0,1");
    }
}
