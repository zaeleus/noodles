//! VCF record genotype field value.

mod array;
pub mod genotype;

pub use self::{array::Array, genotype::Genotype};

use std::str;

/// A VCF record genotype field value.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    /// A 32-bit integer.
    Integer(i32),
    /// A single-precision floating-point.
    Float(f32),
    /// A character.
    Character(char),
    /// A string.
    String(String),
    /// A genotype.
    Genotype(Genotype),
    /// An array.
    Array(Array),
}

impl From<i32> for Value {
    fn from(n: i32) -> Self {
        Self::Integer(n)
    }
}

impl From<f32> for Value {
    fn from(n: f32) -> Self {
        Self::Float(n)
    }
}

impl From<char> for Value {
    fn from(c: char) -> Self {
        Self::Character(c)
    }
}

impl From<&str> for Value {
    fn from(s: &str) -> Self {
        Self::String(s.into())
    }
}

impl From<String> for Value {
    fn from(s: String) -> Self {
        Self::String(s)
    }
}

impl From<Genotype> for Value {
    fn from(genotype: Genotype) -> Self {
        Self::Genotype(genotype)
    }
}

impl From<Vec<Option<i32>>> for Value {
    fn from(values: Vec<Option<i32>>) -> Self {
        Self::Array(Array::Integer(values))
    }
}

impl From<Vec<Option<f32>>> for Value {
    fn from(values: Vec<Option<f32>>) -> Self {
        Self::Array(Array::Float(values))
    }
}

impl From<Vec<Option<char>>> for Value {
    fn from(values: Vec<Option<char>>) -> Self {
        Self::Array(Array::Character(values))
    }
}

impl From<Vec<Option<String>>> for Value {
    fn from(values: Vec<Option<String>>) -> Self {
        Self::Array(Array::String(values))
    }
}

impl<'a> From<&'a Value> for crate::variant::record::samples::series::Value<'a> {
    fn from(value_buf: &'a Value) -> Self {
        match value_buf {
            Value::Integer(n) => Self::Integer(*n),
            Value::Float(n) => Self::Float(*n),
            Value::Character(c) => Self::Character(*c),
            Value::String(s) => Self::String(s.as_ref()),
            Value::Genotype(genotype) => Self::Genotype(Box::new(genotype)),
            Value::Array(array) => Self::Array(array.into()),
        }
    }
}
