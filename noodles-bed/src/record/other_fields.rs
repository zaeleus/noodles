use std::{fmt, iter};

use bstr::{BStr, ByteSlice};

use super::Fields;
use crate::feature::record::other_fields::Value;

/// BED record other fields.
#[derive(Eq, PartialEq)]
pub struct OtherFields<'a, const N: usize>(&'a Fields<N>);

impl<'a, const N: usize> OtherFields<'a, N> {
    pub(super) fn new(fields: &'a Fields<N>) -> Self {
        Self(fields)
    }

    /// Returns whether there are any other fields.
    pub fn is_empty(&self) -> bool {
        self.0.bounds.other_fields_ends.is_empty()
    }

    /// Return the number of other fields.
    pub fn len(&self) -> usize {
        self.0.bounds.other_fields_ends.len()
    }

    /// Returns an other field at the given index.
    pub fn get(&self, i: usize) -> Option<&[u8]> {
        self.0.get(i)
    }

    /// Returns an iterator over other fields.
    pub fn iter(&self) -> impl Iterator<Item = &BStr> {
        let mut i = 0;

        iter::from_fn(move || {
            let field = self.get(i).map(|buf| buf.as_bstr())?;
            i += 1;
            Some(field)
        })
    }
}

impl<'a, const N: usize> fmt::Debug for OtherFields<'a, N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl<'a, const N: usize> crate::feature::record::OtherFields for OtherFields<'a, N> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = Value<'_>> + '_> {
        Box::new(self.iter().map(Value::String))
    }
}
