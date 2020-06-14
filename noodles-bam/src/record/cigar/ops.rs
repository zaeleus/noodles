use std::{convert::TryFrom, mem};

use super::Op;

/// An iterator over the operations of a CIGAR.
///
/// This is created by calling [`Cigar::ops`].
///
/// [`Cigar::ops`]: struct.Cigar.html#method.ops
pub struct Ops<'a> {
    cigar: &'a [u8],
    i: usize,
}

impl<'a> Ops<'a> {
    pub(crate) fn new(cigar: &'a [u8]) -> Self {
        Self { cigar, i: 0 }
    }
}

impl<'a> Iterator for Ops<'a> {
    type Item = Op;

    fn next(&mut self) -> Option<Self::Item> {
        let size = mem::size_of::<u32>();
        let start = self.i * size;

        if start < self.cigar.len() {
            let end = start + size;

            let data = &self.cigar[start..end];
            let op = Op::try_from(data).unwrap();

            self.i += 1;

            Some(op)
        } else {
            None
        }
    }
}
