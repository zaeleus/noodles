use std::io;

use noodles_vcf as vcf;

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
}

impl AsRef<[u8]> for AlternateBases<'_> {
    fn as_ref(&self) -> &[u8] {
        self.src
    }
}

impl vcf::variant::record::AlternateBases for AlternateBases<'_> {
    fn is_empty(&self) -> bool {
        self.len == 0
    }

    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<&str>> + '_> {
        use super::{Value, value::read_value};

        let mut src = self.src;

        Box::new((0..self.len()).map(move |_| match read_value(&mut src)? {
            Some(Value::String(Some(value))) => Ok(value),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid alt value",
            )),
        }))
    }
}
