use std::{iter, mem};

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
