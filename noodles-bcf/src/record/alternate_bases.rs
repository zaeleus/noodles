use std::io;

/// BCF record alternate bases.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlternateBases<'a> {
    src: &'a [u8],
    len: usize,
}

impl<'a> AlternateBases<'a> {
    pub(super) fn new(src: &'a [u8], len: usize) -> Self {
        Self { src, len }
    }

    /// Returns whether there are any alternate alleles.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Returns the number of alternate alleles.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns an iterator over alternate alleles.
    pub fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<&'_ str>>> + '_> {
        use super::{value::read_value, Value};

        let mut src = self.src;

        Box::new((0..self.len()).map(move |_| match read_value(&mut src)? {
            Some(Value::String(value)) => Ok(value),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid alt value",
            )),
        }))
    }
}

impl<'a> AsRef<[u8]> for AlternateBases<'a> {
    fn as_ref(&self) -> &[u8] {
        self.src
    }
}
