use std::{iter, mem};

use noodles_gff as gff;

#[derive(Debug, Eq, PartialEq)]
pub enum Value<'r> {
    String(&'r str),
    Array(Vec<&'r str>),
}

impl<'r> Value<'r> {
    /// An iterator over values.
    pub fn iter(&self) -> Box<dyn Iterator<Item = &'r str> + '_> {
        match self {
            Self::String(value) => Box::new(iter::once(*value)),
            Self::Array(values) => Box::new(values.iter().copied()),
        }
    }

    pub(crate) fn push(&mut self, s: &'r str) {
        match self {
            Self::String(t) => {
                let values = vec![t, s];
                mem::swap(self, &mut Self::Array(values));
            }
            Self::Array(array) => array.push(s),
        }
    }
}

impl From<&Value<'_>> for gff::feature::record_buf::attributes::field::Value {
    fn from(value: &Value<'_>) -> Self {
        match value {
            Value::String(value) => Self::from(*value),
            Value::Array(values) => {
                let vs: Vec<_> = values.iter().map(|s| String::from(*s)).collect();
                Self::from(vs)
            }
        }
    }
}
