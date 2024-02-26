mod field;

use std::io;

use noodles_vcf::{self as vcf, header::StringMaps, variant::record::info::field::Value};

use self::field::read_field;

/// BCF record info.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Info<'a> {
    src: &'a [u8],
    field_count: usize,
}

impl<'a> Info<'a> {
    pub(super) fn new(src: &'a [u8], field_count: usize) -> Self {
        Self { src, field_count }
    }

    /// Returns whether there are any info fields.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the number of info fields.
    pub fn len(&self) -> usize {
        self.field_count
    }

    /// Returns the value with the given key.
    pub fn get<'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
        string_maps: &'h StringMaps,
        key: &str,
    ) -> Option<io::Result<Option<Value<'a>>>> {
        for result in self.iter(header, string_maps) {
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

    /// Returns an iterator over info fields.
    pub fn iter<'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
        string_maps: &'h StringMaps,
    ) -> impl Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a {
        let mut src = self.as_ref();

        (0..self.len()).map(move |_| {
            read_field(&mut src, header, string_maps.strings())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }
}

impl<'a> AsRef<[u8]> for Info<'a> {
    fn as_ref(&self) -> &[u8] {
        self.src
    }
}
