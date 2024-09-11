//! Variant record samples value.

pub mod array;
pub mod genotype;

use std::{borrow::Cow, io};

pub use self::{array::Array, genotype::Genotype};

/// A variant record samples value.
#[derive(Debug)]
pub enum Value<'a> {
    /// A 32-bit integer.
    Integer(i32),
    /// A single-precision floating-point.
    Float(f32),
    /// A character.
    Character(char),
    /// A string.
    String(Cow<'a, str>),
    /// A genotype.
    Genotype(Box<dyn Genotype + 'a>),
    /// An array.
    Array(Array<'a>),
}

impl<'a> TryFrom<Value<'a>> for crate::variant::record_buf::samples::sample::Value {
    type Error = io::Error;

    fn try_from(value: Value<'a>) -> Result<Self, Self::Error> {
        match value {
            Value::Integer(n) => Ok(Self::Integer(n)),
            Value::Float(n) => Ok(Self::Float(n)),
            Value::Character(c) => Ok(Self::Character(c)),
            Value::String(s) => Ok(Self::String(s.into())),
            Value::Genotype(genotype) => genotype.as_ref().try_into().map(Self::Genotype),
            Value::Array(array) => array.try_into().map(Self::Array),
        }
    }
}
