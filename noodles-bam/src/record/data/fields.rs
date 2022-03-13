use std::io;

use crate::{reader::record::data::get_field, record::data::Field};

/// An iterator over BAM record data fields.
///
/// This is created by calling [`crate::record::Data::fields`].
pub struct Fields<'a> {
    buf: &'a [u8],
}

impl<'a> Fields<'a> {
    pub(crate) fn new(buf: &'a [u8]) -> Self {
        Self { buf }
    }
}

impl<'a> Iterator for Fields<'a> {
    type Item = io::Result<Field>;

    fn next(&mut self) -> Option<Self::Item> {
        match get_field(&mut self.buf) {
            Ok(Some(field)) => Some(Ok(field)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}
