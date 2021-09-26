use std::io::{self, BufRead};

use crate::{reader::record::data::read_field, record::data::Field};

/// An iterator over BAM record data fields.
///
/// This is created by calling [`crate::record::Data::fields`].
pub struct Fields<R>
where
    R: BufRead,
{
    inner: R,
}

impl<R> Fields<R>
where
    R: BufRead,
{
    pub(crate) fn new(inner: R) -> Self {
        Self { inner }
    }
}

impl<R> Iterator for Fields<R>
where
    R: BufRead,
{
    type Item = io::Result<Field>;

    fn next(&mut self) -> Option<Self::Item> {
        match read_field(&mut self.inner) {
            Ok(Some(field)) => Some(Ok(field)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}
