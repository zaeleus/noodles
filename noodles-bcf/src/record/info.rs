mod field;

use std::io;

use noodles_vcf::{
    self as vcf,
    variant::record::{info::field::Value, Info as _},
};

use self::field::read_field;

/// BCF record info.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Info<'r> {
    src: &'r [u8],
    field_count: usize,
}

impl<'r> Info<'r> {
    pub(super) fn new(src: &'r [u8], field_count: usize) -> Self {
        Self { src, field_count }
    }

    /// Returns the value with the given key.
    pub fn get<'h: 'r>(
        &'r self,
        header: &'h vcf::Header,
        key: &str,
    ) -> Option<io::Result<Option<Value<'r>>>> {
        for result in self.iter(header) {
            match result {
                Ok((k, v)) => {
                    if k == key {
                        return Some(Ok(v));
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        None
    }
}

impl<'r> AsRef<[u8]> for Info<'r> {
    fn as_ref(&self) -> &[u8] {
        self.src
    }
}

impl<'r> vcf::variant::record::Info for Info<'r> {
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn len(&self) -> usize {
        self.field_count
    }

    fn get<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
        key: &str,
    ) -> Option<io::Result<Option<Value<'a>>>> {
        self.get(header, key)
    }

    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
    ) -> Box<dyn Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a> {
        let mut src = self.as_ref();

        Box::new((0..self.len()).map(move |_| {
            read_field(&mut src, header).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        }))
    }
}
