use std::{borrow::Cow, io};

/// A variant record samples field array value.
#[derive(Clone, Debug, PartialEq)]
pub enum Array {
    /// An array of 32-bit integers.
    Integer(Vec<Option<i32>>),
    /// An array of single-precision floating-points.
    Float(Vec<Option<f32>>),
    /// An array of characters.
    Character(Vec<Option<char>>),
    /// An array of strings.
    String(Vec<Option<String>>),
}

impl<'a> From<&'a Array> for crate::variant::record::samples::series::value::Array<'a> {
    fn from(array_buf: &'a Array) -> Self {
        match array_buf {
            Array::Integer(values) => Self::Integer(Box::new(Values::new(values))),
            Array::Float(values) => Self::Float(Box::new(Values::new(values))),
            Array::Character(values) => Self::Character(Box::new(Values::new(values))),
            Array::String(values) => Self::String(Box::new(Values::new(values))),
        }
    }
}

struct Values<'a, N>(&'a [Option<N>]);

impl<'a, N> Values<'a, N> {
    fn new(values: &'a [Option<N>]) -> Self {
        Self(values)
    }
}

impl<'a, N> crate::variant::record::samples::series::value::array::Values<'a, N> for Values<'a, N>
where
    N: Copy,
{
    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<N>>> + '_> {
        Box::new(self.0.iter().copied().map(Ok))
    }
}

impl<'a> crate::variant::record::samples::series::value::array::Values<'a, Cow<'a, str>>
    for Values<'a, String>
{
    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<Cow<'a, str>>>> + '_> {
        Box::new(self.0.iter().map(|s| Ok(s.as_deref().map(Cow::from))))
    }
}
